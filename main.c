#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <stdint.h>
#include <math.h>
#include <fftw3.h>


#define PCM16_NORMALIZATION_FACTOR 32768.0
#define WINDOW_DURATION 0.09 //90 miliseconds 
#define MAGNITUDE_THRESHOLD 0.01  // keep bins >1% of max magnitude in the window


/* ---- Custom functions ---- */
void mix_to_mono(const int16_t *input, int16_t *output, long long frames, int channels)
{
    for (long long i = 0; i < frames; i++) {
        int sum = 0;

        for (int c = 0; c < channels; c++) {
            sum += input[i * channels + c];
        }

        output[i] = sum / channels;
    }
}

void normalize_mono(const int16_t *input, float *output, long long frames){
    for (long long i = 0; i < frames; i++){
        output[i] = (float)input[i] / PCM16_NORMALIZATION_FACTOR;
    }
}

int compute_quantity_of_windows(int n, int w, int h){ //N = total samples, W = window length (samples), H = hop size (samples)
    if (h == 0){
        printf("Division (hop number) by 0");

        return 0;
    }
    if (n > w) {
        return ((n - w) / h) + 1;
    }

    return 0;
}

void hann_function(float *hann_array, int w) { // w = window length (samples)
    for (int i = 0; i < w; i++) {
        hann_array[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (w - 1)));
    }
}

void fill_windows(const float *normalized_mono, float *windows, const float *hann_array, int window_length, 
    int hop_size, int number_of_windows) 
{
    for (int k = 0; k < number_of_windows; k++) {
        int start = k * hop_size;  // starting index in normalized_mono

        for (int n = 0; n < window_length; n++) {
            windows[k * window_length + n] = normalized_mono[start + n] * hann_array[n];
        }
    }
}


void reconstruct_harmonics_to_wav(fftw_complex **spectra,int window_length,int hop_size,int window_quantity,
        int bins_per_window,  int sample_rate, long long total_samples)
{
    // Precompute max magnitude per bin across all windows
    double *max_magnitude_per_bin = calloc(bins_per_window, sizeof(double));
    for (int bin = 0; bin < bins_per_window; bin++) {
        for (int k = 0; k < window_quantity; k++) {
            double mag = sqrt(spectra[k][bin][0]*spectra[k][bin][0] +
                              spectra[k][bin][1]*spectra[k][bin][1]);
            if (mag > max_magnitude_per_bin[bin])
                max_magnitude_per_bin[bin] = mag;
        }
    }

    // Loop over each harmonic/frequency bin
    for (int bin = 0; bin < bins_per_window; bin++) {

        // Skip bins below threshold
        if (max_magnitude_per_bin[bin] < MAGNITUDE_THRESHOLD)
            continue;

        // Allocate output buffer for this harmonic
        float *harmonic_signal = calloc(total_samples, sizeof(float));
        if (!harmonic_signal) {
            fprintf(stderr, "Memory allocation failed for harmonic %d\n", bin);
            continue;
        }

        // Reconstruct the signal window by window
        for (int k = 0; k < window_quantity; k++) {
            double amplitude = sqrt(spectra[k][bin][0]*spectra[k][bin][0] +
                                    spectra[k][bin][1]*spectra[k][bin][1]);
            double phase = atan2(spectra[k][bin][1], spectra[k][bin][0]);

            for (int n = 0; n < window_length; n++) {
                int pos = k * hop_size + n;
                if (pos >= total_samples) break;

                harmonic_signal[pos] += (float)(amplitude * sin(2.0 * M_PI * bin * n / window_length + phase));
            }
        }

        // Normalize the harmonic signal to avoid clipping
        float max_val = 0.0f;
        for (long long i = 0; i < total_samples; i++) {
            if (fabs(harmonic_signal[i]) > max_val) max_val = fabs(harmonic_signal[i]);
        }
        if (max_val > 0.0f) {
            for (long long i = 0; i < total_samples; i++) {
                harmonic_signal[i] /= max_val;
            }
        }

        // Prepare WAV file info
        SF_INFO sfinfo;
        sfinfo.samplerate = sample_rate;
        sfinfo.channels = 1;
        sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;

        char filename[64];
        snprintf(filename, sizeof(filename), "./audio/deconstructed/harmonic_%d.wav", bin);

        SNDFILE *out_file = sf_open(filename, SFM_WRITE, &sfinfo);
        if (!out_file) {
            fprintf(stderr, "Error opening output file %s: %s\n", filename, sf_strerror(NULL));
            free(harmonic_signal);
            continue;
        }

        // Convert float [-1,1] to PCM16 and write
        int16_t *pcm_buffer = malloc(sizeof(int16_t) * total_samples);
        for (long long i = 0; i < total_samples; i++) {
            pcm_buffer[i] = (int16_t)(harmonic_signal[i] * 32767.0f);
        }

        sf_writef_short(out_file, pcm_buffer, total_samples);
        sf_close(out_file);

        free(pcm_buffer);
        free(harmonic_signal);
    }

    free(max_magnitude_per_bin);

    printf("Significant harmonics written to WAV files.\n");
}



int main(int argc, char *argv[])
{
    if (argc < 2) {
        printf("Usage: %s <audio_file.wav>\n", argv[0]);
        return 1;
    }

    const char *filepath = argv[1];

    SF_INFO info;
    info.format = 0;  /* important: must be zeroed before opening */

    SNDFILE *file = sf_open(filepath, SFM_READ, &info);
    if (file == NULL) {
        printf("Error opening file: %s\n", sf_strerror(NULL));
        return 1;
    }

    /* ---- Validation ---- */

    if ((info.format & SF_FORMAT_TYPEMASK) != SF_FORMAT_WAV) {
        printf("Not a WAV file.\n");
        sf_close(file);
        return 1;
    }

    if ((info.format & SF_FORMAT_SUBMASK) != SF_FORMAT_PCM_16 &&
        (info.format & SF_FORMAT_SUBMASK) != SF_FORMAT_PCM_24 &&
        (info.format & SF_FORMAT_SUBMASK) != SF_FORMAT_PCM_32) {
        printf("WAV file is not linear PCM.\n");
        sf_close(file);
        return 1;
    }

    int16_t *buffer = malloc(info.frames * info.channels * sizeof(int16_t));
    int16_t *mono   = malloc(info.frames * sizeof(int16_t));

    if (!buffer || !mono) {
    printf("Memory allocation failed.\n");
    sf_close(file);
    return 1;
    }

    sf_count_t frames_read = sf_readf_short(file, buffer, info.frames);
    if (frames_read != info.frames) {
        printf("Warning: not all frames were read.\n");
    }

    mix_to_mono(buffer, mono, frames_read, info.channels);


    float *normalized_mono = malloc(info.frames * sizeof(float));


    if (!normalized_mono) {
    printf("Memory allocation failed.\n");
    sf_close(file);
    free(buffer);
    free(mono);
    return 1;
    }

    normalize_mono(mono, normalized_mono, frames_read);

    int window_length = round(info.samplerate * WINDOW_DURATION);
    int hop_size = window_length / 4;

    int window_quantity = compute_quantity_of_windows(frames_read, window_length, hop_size);


    printf("Window length          : %d\n", window_length);
    printf("Hop size               : %d\n", hop_size);
    printf("Window Quantity        : %d\n", window_quantity);

    float *windows = malloc(window_quantity * window_length * sizeof(float));
    float *hann_array = malloc(window_length * sizeof(float));

    if (!windows || !hann_array) {
    printf("Memory allocation failed.\n");
    sf_close(file);
    free(buffer);
    free(mono);
    free(normalized_mono);
    return 1;
    }

    hann_function(hann_array, window_length);
    fill_windows(normalized_mono, windows, hann_array, window_length, hop_size, window_quantity);

    
    int bins_per_window = window_length / 2 + 1;

    fftw_complex **spectra = malloc(sizeof(fftw_complex*) * window_quantity);
    if (!spectra) {
        printf("Memory allocation failed for spectra.\n");
        return 1;
    }
    for (int k = 0; k < window_quantity; k++) {
        spectra[k] = malloc(sizeof(fftw_complex) * bins_per_window);
        if (!spectra[k]) {
            printf("Memory allocation failed for spectra[%d].\n", k);
            return 1;
        }
    }

    double *fft_in = fftw_malloc(sizeof(double) * window_length);
    fftw_complex *fft_out = fftw_malloc(sizeof(fftw_complex) * bins_per_window);

    fftw_plan plan = fftw_plan_dft_r2c_1d(window_length, fft_in, fft_out, FFTW_MEASURE);

    for (int k = 0; k < window_quantity; k++) {

        float *current_window = &windows[k * window_length];

        for (int n = 0; n < window_length; n++) {
            fft_in[n] = (double)current_window[n];
        }

        fftw_execute(plan);

        for (int bin = 0; bin < bins_per_window; bin++) {
            spectra[k][bin][0] = fft_out[bin][0]; // real
            spectra[k][bin][1] = fft_out[bin][1]; // imag
        }
    }

    // spectra[k][bin] now contains all FFT results

    reconstruct_harmonics_to_wav(spectra,window_length,hop_size,window_quantity,bins_per_window,info.samplerate,frames_read);



    /* ---- Metadata output ---- */

    printf("File opened successfully.\n");
    printf("Sample rate       : %d Hz\n", info.samplerate);
    printf("Channels          : %d\n", info.channels);
    printf("Frames            : %lld\n", (long long)info.frames);
    printf("Duration          : %.2f seconds\n",
           (double)info.frames / info.samplerate);

    printf("Format (raw)      : 0x%x\n", info.format);

    sf_close(file);
    free(buffer);
    free(mono);
    free(normalized_mono);
    free(windows);
    free(hann_array);
    fftw_destroy_plan(plan);
    fftw_free(fft_in);
    fftw_free(fft_out);
    for (int k = 0; k < window_quantity; k++) {
        free(spectra[k]);
    }
    free(spectra);
    return 0;
}

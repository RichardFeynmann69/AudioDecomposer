#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <stdint.h>
#include <math.h>

#define PCM16_NORMALIZATION_FACTOR 32768.0
#define WINDOW_DURATION 0.09 //90 miliseconds 


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
    int hop_size = window_length / 2;

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

    printf("\nFirst 10 windows (first 10 samples each):\n");

int max_windows_to_print = window_quantity < 10 ? window_quantity : 10;
int max_samples_to_print = window_length < 10 ? window_length : 10;

for (int k = 0; k < max_windows_to_print; k++) {
    printf("Window %d: ", k);
    for (int n = 0; n < max_samples_to_print; n++) {
        printf("%f ", windows[k * window_length + n]);
    }
    printf("\n");
}



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
    return 0;
}

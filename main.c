#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <stdint.h>

#define PCM16_NORMALIZATION_FACTOR 32768.0


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

    printf("First 10 normalized mono samples:\n");
    for (int i = 0; i < 10 && i < frames_read; i++) {
        printf("%f\n", normalized_mono[i]);
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
    return 0;
}

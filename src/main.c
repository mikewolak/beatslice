/*
 * main.c - Beatslice Command Line Tool
 * 
 * A modern command-line audio beat detection and slicing tool
 * Uses libsndfile for audio I/O and FFTW3 for spectral analysis
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <libgen.h>
#include "beatslice.h"
#include "audio_analysis.h"
#include "fileio.h"
#include "utils.h"

/* Default configuration values */
#define DEFAULT_SENSITIVITY     0.5
#define DEFAULT_MIN_BPM        60.0
#define DEFAULT_MAX_BPM        200.0
#define DEFAULT_FFT_SIZE       2048
#define DEFAULT_HOP_SIZE       512
#define DEFAULT_OUTPUT_FORMAT  "wav"

/* Global configuration structure */
static BeatsliceConfig config = {
    .input_file = NULL,
    .output_prefix = NULL,
    .output_format = DEFAULT_OUTPUT_FORMAT,
    .sensitivity = DEFAULT_SENSITIVITY,
    .min_bpm = DEFAULT_MIN_BPM,
    .max_bpm = DEFAULT_MAX_BPM,
    .fft_size = DEFAULT_FFT_SIZE,
    .hop_size = DEFAULT_HOP_SIZE,
    .use_adaptive_threshold = 1,
    .use_zero_crossings = 0,
    .verbose = 0,
    .export_analysis = 0,
    .min_slice_length = 0.1,  /* 100ms minimum slice */
    .normalize_output = 0
};

void print_usage(const char *program_name) 
{
    printf("Beatslice - Audio Beat Detection and Slicing Tool\n");
    printf("Usage: %s [OPTIONS] <input_file>\n\n", program_name);
    
    printf("Required:\n");
    printf("  <input_file>              Input audio file (wav, aiff, flac, mp3, etc.)\n\n");
    
    printf("Options:\n");
    printf("  -o, --output PREFIX       Output file prefix (default: input filename)\n");
    printf("  -f, --format FORMAT       Output format: wav, aiff, flac (default: wav)\n");
    printf("  -s, --sensitivity FLOAT   Detection sensitivity 0.0-1.0 (default: 0.5)\n");
    printf("  -b, --min-bpm FLOAT       Minimum BPM for detection (default: 60)\n");
    printf("  -B, --max-bpm FLOAT       Maximum BPM for detection (default: 200)\n");
    printf("  -w, --window-size INT     FFT window size, power of 2 (default: 2048)\n");
    printf("  -h, --hop-size INT        FFT hop size (default: 512)\n");
    printf("  -a, --adaptive            Use adaptive thresholding (default: on)\n");
    printf("  -z, --zero-crossings      Include zero-crossing analysis\n");
    printf("  -m, --min-length FLOAT    Minimum slice length in seconds (default: 0.1)\n");
    printf("  -n, --normalize           Normalize output slices\n");
    printf("  -e, --export-analysis     Export analysis data to JSON\n");
    printf("  -v, --verbose             Verbose output\n");
    printf("  --help                    Show this help message\n\n");
    
    printf("Examples:\n");
    printf("  %s song.wav\n", program_name);
    printf("  %s -s 0.8 -o drums drum_loop.flac\n", program_name);
    printf("  %s --sensitivity 0.3 --format aiff --verbose track.wav\n", program_name);
    printf("  %s -s 0.7 -b 120 -B 140 --export-analysis dubstep.wav\n", program_name);
}

void print_version(void)
{
    printf("Beatslice v1.0.0\n");
    printf("Built with libsndfile %s, FFTW3\n", sf_version_string());
}

int parse_arguments(int argc, char **argv)
{
    int opt;
    int option_index = 0;
    
    static struct option long_options[] = {
        {"output",          required_argument, 0, 'o'},
        {"format",          required_argument, 0, 'f'},
        {"sensitivity",     required_argument, 0, 's'},
        {"min-bpm",         required_argument, 0, 'b'},
        {"max-bpm",         required_argument, 0, 'B'},
        {"window-size",     required_argument, 0, 'w'},
        {"hop-size",        required_argument, 0, 'h'},
        {"adaptive",        no_argument,       0, 'a'},
        {"zero-crossings",  no_argument,       0, 'z'},
        {"min-length",      required_argument, 0, 'm'},
        {"normalize",       no_argument,       0, 'n'},
        {"export-analysis", no_argument,       0, 'e'},
        {"verbose",         no_argument,       0, 'v'},
        {"help",            no_argument,       0, '?'},
        {"version",         no_argument,       0, 'V'},
        {0, 0, 0, 0}
    };
    
    while ((opt = getopt_long(argc, argv, "o:f:s:b:B:w:h:azm:nevV?", 
                              long_options, &option_index)) != -1) {
        switch (opt) {
            case 'o':
                config.output_prefix = strdup(optarg);
                break;
                
            case 'f':
                config.output_format = strdup(optarg);
                if (strcmp(config.output_format, "wav") != 0 && 
                    strcmp(config.output_format, "aiff") != 0 && 
                    strcmp(config.output_format, "flac") != 0) {
                    fprintf(stderr, "Error: Unsupported output format '%s'\n", 
                           config.output_format);
                    return -1;
                }
                break;
                
            case 's':
                config.sensitivity = atof(optarg);
                if (config.sensitivity < 0.0 || config.sensitivity > 1.0) {
                    fprintf(stderr, "Error: Sensitivity must be between 0.0 and 1.0\n");
                    return -1;
                }
                break;
                
            case 'b':
                config.min_bpm = atof(optarg);
                if (config.min_bpm <= 0) {
                    fprintf(stderr, "Error: Minimum BPM must be positive\n");
                    return -1;
                }
                break;
                
            case 'B':
                config.max_bpm = atof(optarg);
                if (config.max_bpm <= 0) {
                    fprintf(stderr, "Error: Maximum BPM must be positive\n");
                    return -1;
                }
                break;
                
            case 'w':
                config.fft_size = atoi(optarg);
                if (!is_power_of_two(config.fft_size) || config.fft_size < 256) {
                    fprintf(stderr, "Error: FFT size must be a power of 2, >= 256\n");
                    return -1;
                }
                break;
                
            case 'h':
                config.hop_size = atoi(optarg);
                if (config.hop_size <= 0) {
                    fprintf(stderr, "Error: Hop size must be positive\n");
                    return -1;
                }
                break;
                
            case 'a':
                config.use_adaptive_threshold = 1;
                break;
                
            case 'z':
                config.use_zero_crossings = 1;
                break;
                
            case 'm':
                config.min_slice_length = atof(optarg);
                if (config.min_slice_length <= 0) {
                    fprintf(stderr, "Error: Minimum slice length must be positive\n");
                    return -1;
                }
                break;
                
            case 'n':
                config.normalize_output = 1;
                break;
                
            case 'e':
                config.export_analysis = 1;
                break;
                
            case 'v':
                config.verbose = 1;
                break;
                
            case 'V':
                print_version();
                exit(0);
                break;
                
            case '?':
            default:
                print_usage(argv[0]);
                exit(0);
                break;
        }
    }
    
    /* Check for required input file */
    if (optind >= argc) {
        fprintf(stderr, "Error: No input file specified\n");
        print_usage(argv[0]);
        return -1;
    }
    
    config.input_file = strdup(argv[optind]);
    
    /* Set default output prefix if not specified */
    if (!config.output_prefix) {
        char *input_copy = strdup(config.input_file);
        char *base_name = basename(input_copy);
        char *dot = strrchr(base_name, '.');
        if (dot) *dot = '\0';  /* Remove extension */
        config.output_prefix = strdup(base_name);
        free(input_copy);
    }
    
    /* Validate BPM range */
    if (config.min_bpm >= config.max_bpm) {
        fprintf(stderr, "Error: Minimum BPM must be less than maximum BPM\n");
        return -1;
    }
    
    return 0;
}

void print_config(void)
{
    if (!config.verbose) return;
    
    printf("Configuration:\n");
    printf("  Input file: %s\n", config.input_file);
    printf("  Output prefix: %s\n", config.output_prefix);
    printf("  Output format: %s\n", config.output_format);
    printf("  Sensitivity: %.2f\n", config.sensitivity);
    printf("  BPM range: %.1f - %.1f\n", config.min_bpm, config.max_bpm);
    printf("  FFT size: %d, Hop size: %d\n", config.fft_size, config.hop_size);
    printf("  Adaptive threshold: %s\n", config.use_adaptive_threshold ? "yes" : "no");
    printf("  Zero crossings: %s\n", config.use_zero_crossings ? "yes" : "no");
    printf("  Min slice length: %.2fs\n", config.min_slice_length);
    printf("  Normalize output: %s\n", config.normalize_output ? "yes" : "no");
    printf("\n");
}

int main(int argc, char **argv)
{
    AudioData audio_data;
    AnalysisResults analysis;
    int result = 0;
    
    /* Parse command line arguments */
    if (parse_arguments(argc, argv) != 0) {
        return 1;
    }
    
    print_config();
    
    /* Initialize audio analysis system */
    if (analysis_init(&config) != 0) {
        fprintf(stderr, "Error: Failed to initialize analysis system\n");
        return 1;
    }
    
    /* Load input audio file */
    if (config.verbose) {
        printf("Loading audio file: %s\n", config.input_file);
    }
    
    if (load_audio_file(config.input_file, &audio_data) != 0) {
        fprintf(stderr, "Error: Failed to load audio file\n");
        result = 1;
        goto cleanup;
    }
    
    if (config.verbose) {
        printf("Audio loaded: %.2fs, %d Hz, %d channels\n", 
               audio_data.duration, audio_data.sample_rate, audio_data.channels);
    }
    
    /* Perform beat detection analysis */
    if (config.verbose) {
        printf("Performing beat analysis...\n");
    }
    
    if (perform_beat_analysis(&audio_data, &config, &analysis) != 0) {
        fprintf(stderr, "Error: Beat analysis failed\n");
        result = 1;
        goto cleanup;
    }
    
    if (config.verbose) {
        printf("Analysis complete: %d beats detected, estimated BPM: %.1f\n", 
               analysis.beat_count, analysis.estimated_bpm);
    }
    
    /* Export beat slices */
    if (config.verbose) {
        printf("Exporting beat slices...\n");
    }
    
    int slices_exported = export_beat_slices(&audio_data, &analysis, &config);
    if (slices_exported <= 0) {
        fprintf(stderr, "Error: Failed to export slices\n");
        result = 1;
        goto cleanup;
    }
    
    printf("Successfully exported %d slices\n", slices_exported);
    
    /* Export analysis data if requested */
    if (config.export_analysis) {
        if (config.verbose) {
            printf("Exporting analysis data...\n");
        }
        
        if (export_analysis_data(&analysis, &config) != 0) {
            fprintf(stderr, "Warning: Failed to export analysis data\n");
        }
    }

cleanup:
    /* Clean up resources */
    if (audio_data.samples) {
        free(audio_data.samples);
    }
    
    if (analysis.beat_times) {
        free(analysis.beat_times);
    }
    
    if (analysis.confidence) {
        free(analysis.confidence);
    }
    
    if (analysis.energy_envelope) {
        free(analysis.energy_envelope);
    }
    
    analysis_cleanup();
    
    if (config.input_file) free((char*)config.input_file);
    if (config.output_prefix) free((char*)config.output_prefix);
    
    return result;
}

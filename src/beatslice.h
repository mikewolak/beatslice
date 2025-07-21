/*
 * beatslice.h - Main header for Beatslice audio analysis tool
 * 
 * Defines core data structures and function prototypes for beat detection
 * and audio slicing using libsndfile and FFTW3.
 */

#ifndef BEATSLICE_H
#define BEATSLICE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sndfile.h>
#include <fftw3.h>

/* Version information */
#define BEATSLICE_VERSION_MAJOR  1
#define BEATSLICE_VERSION_MINOR  0
#define BEATSLICE_VERSION_PATCH  0

/* Analysis constants */
#define MAX_BEATS               2048
#define MAX_ANALYSIS_FRAMES     16384
#define MIN_SLICE_DURATION      0.01   /* 10ms minimum slice */
#define MAX_SLICE_DURATION      30.0   /* 30s maximum slice */

/* Configuration structure */
typedef struct {
    const char *input_file;
    const char *output_prefix;
    const char *output_format;
    double sensitivity;              /* 0.0 - 1.0 */
    double min_bpm;                 /* Minimum BPM constraint */
    double max_bpm;                 /* Maximum BPM constraint */
    int fft_size;                   /* FFT window size */
    int hop_size;                   /* FFT hop size */
    int use_adaptive_threshold;     /* Boolean flag */
    int use_zero_crossings;         /* Boolean flag */
    int verbose;                    /* Verbose output */
    int export_analysis;            /* Export JSON analysis data */
    double min_slice_length;        /* Minimum slice length in seconds */
    int normalize_output;           /* Normalize output slices */
} BeatsliceConfig;

/* Audio data structure */
typedef struct {
    float *samples;                 /* Interleaved audio samples */
    sf_count_t frame_count;         /* Number of frames */
    int sample_rate;                /* Sample rate in Hz */
    int channels;                   /* Number of channels */
    double duration;                /* Duration in seconds */
    SF_INFO sf_info;               /* libsndfile info structure */
} AudioData;

/* Spectral analysis frame */
typedef struct {
    double time_position;           /* Time position in seconds */
    double rms_energy;             /* RMS energy for this frame */
    double spectral_centroid;     /* Spectral centroid (brightness) */
    double spectral_rolloff;      /* Spectral rolloff frequency */
    double zero_crossing_rate;    /* Zero crossing rate */
    double low_energy;            /* Energy in low freq band (20-250 Hz) */
    double mid_energy;            /* Energy in mid freq band (250-2000 Hz) */
    double high_energy;           /* Energy in high freq band (2000+ Hz) */
    double onset_strength;        /* Combined onset detection function */
} SpectralFrame;

/* Beat detection results */
typedef struct {
    double *beat_times;            /* Array of beat times in seconds */
    double *confidence;            /* Confidence scores for each beat */
    int beat_count;               /* Number of beats detected */
    double estimated_bpm;         /* Estimated BPM */
    double tempo_confidence;      /* Confidence in BPM estimation */
    
    /* Analysis data for export */
    SpectralFrame *frames;        /* Spectral analysis frames */
    int frame_count;             /* Number of analysis frames */
    double *energy_envelope;     /* Low-frequency energy envelope */
    double *onset_function;      /* Onset detection function */
    double analysis_hop_time;    /* Time between analysis frames */
} AnalysisResults;

/* Enhanced sensitivity parameters */
typedef struct {
    double peak_threshold;        /* Base peak detection threshold */
    double adaptive_factor;       /* Adaptive threshold scaling */
    double min_peak_distance;    /* Minimum time between peaks */
    double onset_sensitivity;    /* Onset detection sensitivity */
    double spectral_flux_weight; /* Weight for spectral flux */
    double energy_weight;        /* Weight for energy changes */
    double hfc_weight;           /* Weight for high frequency content */
} SensitivityParams;

/* Beat detection configuration */
typedef struct {
    SensitivityParams sensitivity;
    double min_tempo;            /* Minimum tempo in BPM */
    double max_tempo;            /* Maximum tempo in BPM */
    int use_adaptive_threshold;  /* Use adaptive thresholding */
    int use_multiple_features;   /* Combine multiple onset features */
    int median_filter_length;    /* Median filter for onset function */
    double onset_combine_time;   /* Time window for combining onsets */
} BeatDetectionConfig;

/* Slice export information */
typedef struct {
    double start_time;           /* Start time in seconds */
    double end_time;            /* End time in seconds */
    sf_count_t start_frame;     /* Start frame index */
    sf_count_t frame_count;     /* Number of frames in slice */
    double confidence;          /* Beat detection confidence */
    char filename[512];         /* Output filename */
} SliceInfo;

/*----- Function Prototypes -----*/

/* Main analysis functions */
int analysis_init(const BeatsliceConfig *config);
void analysis_cleanup(void);
int perform_beat_analysis(const AudioData *audio, const BeatsliceConfig *config, 
                         AnalysisResults *results);

/* Audio I/O functions (fileio.c) */
int load_audio_file(const char *filename, AudioData *audio_data);
int export_beat_slices(const AudioData *audio, const AnalysisResults *analysis, 
                      const BeatsliceConfig *config);
int export_analysis_data(const AnalysisResults *analysis, const BeatsliceConfig *config);
const char* get_output_filename(const char *prefix, int slice_index, 
                               const char *format);

/* Enhanced audio analysis functions (audio_analysis.c) */
int compute_spectral_features(const AudioData *audio, const BeatsliceConfig *config, 
                             SpectralFrame *frames, int *frame_count);

/* Error handling */
typedef enum {
    BEATSLICE_SUCCESS = 0,
    BEATSLICE_ERROR_FILE_NOT_FOUND,
    BEATSLICE_ERROR_UNSUPPORTED_FORMAT,
    BEATSLICE_ERROR_MEMORY_ALLOCATION,
    BEATSLICE_ERROR_INVALID_PARAMETER,
    BEATSLICE_ERROR_ANALYSIS_FAILED,
    BEATSLICE_ERROR_EXPORT_FAILED
} BeatsliceError;

const char* beatslice_error_string(BeatsliceError error);

/* Debug and logging */
#ifdef DEBUG
#define LOG_DEBUG(fmt, ...) fprintf(stderr, "[DEBUG] " fmt "\n", ##__VA_ARGS__)
#else
#define LOG_DEBUG(fmt, ...)
#endif

#define LOG_INFO(fmt, ...) do { if (1) fprintf(stdout, "[INFO] " fmt "\n", ##__VA_ARGS__); } while(0)
#define LOG_WARNING(fmt, ...) fprintf(stderr, "[WARNING] " fmt "\n", ##__VA_ARGS__)
#define LOG_ERROR(fmt, ...) fprintf(stderr, "[ERROR] " fmt "\n", ##__VA_ARGS__)

#endif /* BEATSLICE_H */

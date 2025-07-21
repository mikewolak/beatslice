/*
 * audio_analysis.h - Enhanced Audio Analysis Header
 * 
 * Function declarations and data structures for advanced beat detection,
 * spectral analysis, and onset detection using FFTW3.
 */

#ifndef AUDIO_ANALYSIS_H
#define AUDIO_ANALYSIS_H

#include <fftw3.h>
#include "beatslice.h"

/* Analysis algorithm constants */
#define DEFAULT_FFT_SIZE            2048
#define DEFAULT_HOP_SIZE            512
#define DEFAULT_WINDOW_TYPE         WINDOW_HANN
#define DEFAULT_OVERLAP_FACTOR      0.75

/* Spectral analysis parameters */
#define MIN_FFT_SIZE               256
#define MAX_FFT_SIZE               16384
#define MIN_HOP_SIZE               64
#define MAX_HOP_SIZE               4096

/* Onset detection parameters */
#define MIN_ONSET_THRESHOLD        0.01
#define MAX_ONSET_THRESHOLD        1.0
#define DEFAULT_ONSET_THRESHOLD    0.3
#define DEFAULT_MEDIAN_FILTER_SIZE 5
#define MAX_MEDIAN_FILTER_SIZE     15

/* Beat detection constraints */
#define MIN_BEAT_INTERVAL_MS       200    /* 300 BPM maximum */
#define MAX_BEAT_INTERVAL_MS       2000   /* 30 BPM minimum */
#define DEFAULT_MIN_BPM            60.0
#define DEFAULT_MAX_BPM            200.0

/* Frequency band definitions (Hz) */
#define BASS_FREQ_MIN              20.0
#define BASS_FREQ_MAX              250.0
#define MID_FREQ_MIN               250.0
#define MID_FREQ_MAX               2000.0
#define HIGH_FREQ_MIN              2000.0
#define HIGH_FREQ_MAX              8000.0

/* Window function types */
typedef enum {
    WINDOW_HANN = 0,
    WINDOW_HAMMING,
    WINDOW_BLACKMAN,
    WINDOW_RECTANGULAR,
    WINDOW_KAISER
} WindowType;

/* Onset detection methods */
typedef enum {
    ONSET_ENERGY = 0,           /* Energy-based detection */
    ONSET_SPECTRAL_FLUX,        /* Spectral flux method */
    ONSET_HIGH_FREQUENCY,       /* High frequency content */
    ONSET_COMPLEX_DOMAIN,       /* Complex domain method */
    ONSET_COMBINED              /* Combination of methods */
} OnsetMethod;

/* Beat detection algorithms */
typedef enum {
    BEAT_SIMPLE_PEAKS = 0,      /* Simple peak picking */
    BEAT_ADAPTIVE_THRESHOLD,    /* Adaptive threshold method */
    BEAT_DYNAMIC_PROGRAMMING,   /* Dynamic programming approach */
    BEAT_TEMPLATE_MATCHING      /* Template matching method */
} BeatDetectionMethod;

/* Advanced spectral features */
typedef struct {
    double spectral_centroid;      /* Frequency center of mass */
    double spectral_rolloff;       /* 85% energy cutoff frequency */
    double spectral_flatness;     /* Spectral flatness measure */
    double spectral_flux;          /* Change in spectral content */
    double zero_crossing_rate;     /* Zero crossing rate */
    double mfcc[13];              /* Mel-frequency cepstral coefficients */
    double chroma[12];            /* Chromagram features */
    double tonnetz[6];            /* Tonal centroid features */
} AdvancedSpectralFeatures;

/* Enhanced onset detection configuration */
typedef struct {
    OnsetMethod method;            /* Primary onset detection method */
    double threshold;              /* Base detection threshold */
    double threshold_factor;       /* Adaptive threshold scaling */
    int median_filter_length;     /* Median filter size for smoothing */
    int use_adaptive_threshold;   /* Enable adaptive thresholding */
    int combine_methods;          /* Combine multiple onset methods */
    double spectral_flux_weight;  /* Weight for spectral flux */
    double energy_weight;         /* Weight for energy changes */
    double hfc_weight;           /* Weight for high frequency content */
    double whitening_factor;      /* Spectral whitening factor */
} OnsetDetectionConfig;

/* Enhanced beat detection configuration */
typedef struct {
    BeatDetectionMethod method;    /* Beat detection algorithm */
    SensitivityParams sensitivity; /* Sensitivity parameters */
    double min_tempo;             /* Minimum tempo in BPM */
    double max_tempo;             /* Maximum tempo in BPM */
    int use_tempo_constraints;    /* Enforce tempo constraints */
    int use_adaptive_threshold;   /* Use adaptive thresholding */
    int use_multiple_features;    /* Combine multiple onset methods */
    int median_filter_length;     /* Median filter for onset function */
    double onset_combine_time;    /* Time window for combining onsets */
    int peak_interpolation;       /* Use sub-sample peak interpolation */
    double confidence_threshold;  /* Minimum confidence for beat */
} EnhancedBeatDetectionConfig;

/* Analysis state and buffers */
typedef struct {
    /* FFTW3 objects */
    fftw_plan fft_plan;           /* Forward FFT plan */
    double *fft_input;            /* FFT input buffer */
    fftw_complex *fft_output;     /* FFT output buffer */
    
    /* Analysis buffers */
    double *window_function;      /* Window function coefficients */
    double *spectrum_current;     /* Current frame spectrum */
    double *spectrum_previous;    /* Previous frame spectrum */
    double *onset_function;       /* Onset detection function */
    double *adaptive_threshold;   /* Adaptive threshold buffer */
    
    /* Configuration */
    int fft_size;                 /* FFT window size */
    int hop_size;                 /* Hop size in samples */
    int sample_rate;              /* Audio sample rate */
    WindowType window_type;       /* Window function type */
    
    /* State variables */
    int initialized;              /* Whether analysis is initialized */
    int frame_count;             /* Number of processed frames */
    double analysis_time;        /* Current analysis time position */
} AnalysisState;

/*----- Core Analysis Functions -----*/

/**
 * Initialize the audio analysis system
 *
 * Sets up FFTW3 plans, allocates buffers, and configures the analysis
 * system according to the provided configuration.
 *
 * @param config Pointer to beatslice configuration
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int analysis_init(const BeatsliceConfig *config);

/**
 * Clean up audio analysis system
 *
 * Frees all allocated memory, destroys FFTW plans, and resets the
 * analysis system to uninitialized state.
 */
void analysis_cleanup(void);

/**
 * Perform complete beat analysis on audio data
 *
 * Main analysis function that coordinates spectral analysis, onset detection,
 * and beat detection to produce comprehensive analysis results.
 *
 * @param audio Pointer to input audio data
 * @param config Pointer to analysis configuration
 * @param results Pointer to store analysis results
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
/**
 * Perform complete beat analysis on audio data
 *
 * Main analysis function that coordinates spectral analysis, onset detection,
 * and beat detection to produce comprehensive analysis results.
 *
 * @param audio Pointer to input audio data
 * @param config Pointer to analysis configuration
 * @param results Pointer to store analysis results
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int perform_beat_analysis(const AudioData *audio, const BeatsliceConfig *config,
                         AnalysisResults *results);

/**
 * Detect beats using original algorithm
 *
 * Uses the original proven beat detection algorithm adapted for new data structures.
 *
 * @param audio Pointer to input audio data
 * @param config Pointer to analysis configuration
 * @param results Pointer to store analysis results
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int detect_beats_original(const AudioData *audio, const BeatsliceConfig *config,
                         AnalysisResults *results);

/**
 * Compute spectral features for all audio frames
 *
 * Processes the input audio using overlapping FFT windows to extract
 * comprehensive spectral features for each analysis frame.
 *
 * @param audio Pointer to input audio data
 * @param config Pointer to analysis configuration
 * @param frames Array to store computed spectral frames
 * @param frame_count Pointer to store number of frames computed
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int compute_spectral_features(const AudioData *audio, const BeatsliceConfig *config,
                             SpectralFrame *frames, int *frame_count);

/**
 * Compute advanced spectral features
 *
 * Calculates advanced spectral features including MFCCs, chromagram,
 * and tonal centroid features for music information retrieval.
 *
 * @param spectrum Magnitude spectrum of current frame
 * @param spectrum_size Size of spectrum array
 * @param sample_rate Audio sample rate
 * @param features Pointer to store advanced features
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int compute_advanced_spectral_features(const double *spectrum, int spectrum_size,
                                      int sample_rate, AdvancedSpectralFeatures *features);

/**
 * Calculate spectral centroid
 *
 * Computes the spectral centroid (center of mass) of the frequency spectrum.
 *
 * @param spectrum Magnitude spectrum array
 * @param length Length of spectrum array
 * @param sample_rate Audio sample rate
 * @return Spectral centroid in Hz
 */
double calculate_spectral_centroid(const double *spectrum, int length, double sample_rate);

/**
 * Calculate spectral rolloff
 *
 * Computes the spectral rolloff frequency (frequency below which a specified
 * percentage of total energy is contained).
 *
 * @param spectrum Magnitude spectrum array
 * @param length Length of spectrum array
 * @param sample_rate Audio sample rate
 * @param rolloff_point Percentage of energy (e.g., 0.85 for 85%)
 * @return Spectral rolloff frequency in Hz
 */
double calculate_spectral_rolloff(const double *spectrum, int length, 
                                 double sample_rate, double rolloff_point);

/**
 * Calculate spectral flux
 *
 * Computes the spectral flux between two consecutive spectra, measuring
 * the rate of change in spectral content.
 *
 * @param current_spectrum Current frame spectrum
 * @param previous_spectrum Previous frame spectrum
 * @param length Length of spectrum arrays
 * @return Spectral flux value
 */
double calculate_spectral_flux(const double *current_spectrum, 
                              const double *previous_spectrum, int length);

/**
 * Calculate high frequency content
 *
 * Computes a measure of high frequency content by weighting spectrum
 * bins by their frequency index.
 *
 * @param spectrum Magnitude spectrum array
 * @param length Length of spectrum array
 * @return High frequency content value
 */
double calculate_high_frequency_content(const double *spectrum, int length);

/*----- Onset Detection Functions -----*/

/**
 * Detect onsets in spectral frames
 *
 * Processes spectral frames to generate an onset detection function
 * using the specified method and configuration.
 *
 * @param frames Array of spectral frames
 * @param frame_count Number of frames
 * @param config Onset detection configuration
 * @param onset_function Output array for onset detection function
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int detect_onsets(const SpectralFrame *frames, int frame_count,
                 const EnhancedBeatDetectionConfig *config, double *onset_function);

/**
 * Calculate onset strength between two frames
 *
 * Computes a combined onset strength measure using multiple spectral
 * features and weighting factors.
 *
 * @param current Pointer to current spectral frame
 * @param previous Pointer to previous spectral frame
 * @return Combined onset strength value
 */
double calculate_onset_strength(const SpectralFrame *current, 
                               const SpectralFrame *previous);

/**
 * Apply onset detection method
 *
 * Applies a specific onset detection method to compute onset strength
 * from spectral data.
 *
 * @param method Onset detection method to use
 * @param current_spectrum Current frame spectrum
 * @param previous_spectrum Previous frame spectrum
 * @param spectrum_size Size of spectrum arrays
 * @param sample_rate Audio sample rate
 * @return Onset strength value for the method
 */
double apply_onset_method(OnsetMethod method, const double *current_spectrum,
                         const double *previous_spectrum, int spectrum_size,
                         int sample_rate);

/*----- Beat Detection Functions -----*/

/**
 * Enhanced beat detection from onset function
 *
 * Detects beats using advanced algorithms on the onset detection function
 * with tempo constraints and confidence scoring.
 *
 * @param onset_function Onset detection function array
 * @param function_length Length of onset function
 * @param hop_time Time between analysis frames in seconds
 * @param config Beat detection configuration
 * @param results Pointer to store beat detection results
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int detect_beats_enhanced(const double *onset_function, int function_length,
                         double hop_time, const EnhancedBeatDetectionConfig *config,
                         AnalysisResults *results);

/**
 * Find peaks in onset function with advanced constraints
 *
 * Locates peaks in the onset function using adaptive thresholding,
 * temporal constraints, and confidence scoring.
 *
 * @param signal Input signal array
 * @param length Length of signal array
 * @param config Beat detection configuration
 * @param peak_times Output array for peak times
 * @param peak_confidences Output array for peak confidence scores
 * @param max_peaks Maximum number of peaks to detect
 * @param hop_time Time between samples in seconds
 * @return Number of peaks detected
 */
int find_peaks_advanced(const double *signal, int length,
                       const EnhancedBeatDetectionConfig *config,
                       double *peak_times, double *peak_confidences,
                       int max_peaks, double hop_time);

/**
 * Perform complete beat analysis on audio data
 *
 * Main analysis function that coordinates spectral analysis, onset detection,
 * and beat detection to produce comprehensive analysis results.
 *
 * @param audio Pointer to input audio data
 * @param config Pointer to analysis configuration
 * @param results Pointer to store analysis results
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int perform_beat_analysis(const AudioData *audio, const BeatsliceConfig *config,
                         AnalysisResults *results);

/**
 * Detect beats using original algorithm
 *
 * Uses the original proven beat detection algorithm adapted for new data structures.
 *
 * @param audio Pointer to input audio data
 * @param config Pointer to analysis configuration
 * @param results Pointer to store analysis results
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int detect_beats_original(const AudioData *audio, const BeatsliceConfig *config,
                         AnalysisResults *results);

/**
 * Estimate tempo from beat positions
 *
 * Analyzes beat positions to estimate the most likely tempo using
 * histogram analysis and peak detection.
 *
 * @param beat_times Array of beat times in seconds
 * @param beat_count Number of beats
 * @param min_bpm Minimum BPM to consider
 * @param max_bpm Maximum BPM to consider
 * @return Estimated BPM value
 */
double estimate_tempo(const double *beat_times, int beat_count, 
                     double min_bpm, double max_bpm);

/*----- Configuration and Sensitivity Functions -----*/

/**
 * Configure sensitivity parameters
 *
 * Maps user sensitivity setting to detailed algorithm parameters
 * using non-linear response curves and feature weighting.
 *
 * @param sensitivity User sensitivity value (0.0-1.0)
 * @param config Input configuration
 * @param beat_config Output beat detection configuration
 */
void configure_sensitivity(double sensitivity, const BeatsliceConfig *config,
                          EnhancedBeatDetectionConfig *beat_config);

/**
 * Generate adaptive threshold for signal
 *
 * Creates a dynamic threshold that adapts to local signal characteristics
 * using sliding window analysis.
 *
 * @param signal Input signal array
 * @param length Length of signal array
 * @param factor Scaling factor for threshold
 * @param threshold Output threshold array
 */
void adaptive_threshold(const double *signal, int length, double factor,
                       double *threshold);

/**
 * Configure onset detection parameters
 *
 * Sets up onset detection configuration based on audio characteristics
 * and user preferences.
 *
 * @param audio_data Pointer to audio data for analysis
 * @param user_config User configuration
 * @param onset_config Output onset detection configuration
 */
void configure_onset_detection(const AudioData *audio_data,
                              const BeatsliceConfig *user_config,
                              OnsetDetectionConfig *onset_config);

/*----- Window Function Utilities -----*/

/**
 * Generate window function coefficients
 *
 * Generates coefficients for the specified window function type.
 *
 * @param window_type Type of window function
 * @param size Window size in samples
 * @param coefficients Output array for window coefficients
 */
void generate_window_function(WindowType window_type, int size, double *coefficients);

/**
 * Hann window function
 *
 * Computes a single Hann window coefficient.
 *
 * @param n Sample index
 * @param N Total window size
 * @return Window coefficient value
 */
double hann_window(int n, int N);

/**
 * Apply window function to signal
 *
 * Multiplies input signal by window function coefficients.
 *
 * @param signal Input/output signal array
 * @param window Window function coefficients
 * @param length Length of arrays
 */
void apply_window_function(double *signal, const double *window, int length);

/*----- Utility and Helper Functions -----*/

/**
 * Convert frequency to FFT bin index
 *
 * Converts a frequency in Hz to the corresponding FFT bin index.
 *
 * @param frequency Frequency in Hz
 * @param fft_size FFT size
 * @param sample_rate Sample rate in Hz
 * @return FFT bin index
 */
int frequency_to_bin(double frequency, int fft_size, int sample_rate);

/**
 * Convert FFT bin index to frequency
 *
 * Converts an FFT bin index to the corresponding frequency in Hz.
 *
 * @param bin_index FFT bin index
 * @param fft_size FFT size
 * @param sample_rate Sample rate in Hz
 * @return Frequency in Hz
 */
double bin_to_frequency(int bin_index, int fft_size, int sample_rate);

/**
 * Calculate RMS energy of signal segment
 *
 * Computes the Root Mean Square energy of a signal segment.
 *
 * @param samples Input audio samples
 * @param length Number of samples
 * @return RMS energy value
 */
double calculate_rms_energy(const float *samples, int length);

/**
 * Zero-pad signal to power of 2 length
 *
 * Pads a signal with zeros to the next power of 2 length for efficient FFT.
 *
 * @param input Input signal
 * @param input_length Input signal length
 * @param output Output padded signal
 * @param output_length Desired output length (power of 2)
 */
void zero_pad_signal(const double *input, int input_length,
                     double *output, int output_length);

/**
 * Validate analysis parameters
 *
 * Checks that analysis parameters are within valid ranges and
 * compatible with each other.
 *
 * @param config Configuration to validate
 * @return 1 if valid, 0 if invalid
 */
int validate_analysis_parameters(const BeatsliceConfig *config);

/**
 * Get analysis state information
 *
 * Returns pointer to current analysis state for debugging and monitoring.
 *
 * @return Pointer to AnalysisState structure
 */
const AnalysisState* get_analysis_state(void);

/**
 * Reset analysis state
 *
 * Resets internal analysis state while keeping buffers allocated.
 * Useful for processing multiple files without full reinitialization.
 */
void reset_analysis_state(void);

/*----- Debug and Profiling Functions -----*/

/**
 * Export intermediate analysis data
 *
 * Exports intermediate analysis data (spectra, onset function, etc.)
 * for debugging and algorithm development.
 *
 * @param filename Output filename prefix
 * @param analysis_state Current analysis state
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int export_debug_data(const char *filename, const AnalysisState *analysis_state);

/**
 * Print analysis timing information
 *
 * Prints detailed timing information for different analysis stages.
 */
void print_analysis_timing(void);

/**
 * Get analysis performance metrics
 *
 * Returns performance metrics for the analysis process.
 *
 * @param frames_per_second Processing speed in frames per second
 * @param realtime_factor Real-time factor (1.0 = real-time)
 * @param memory_usage Memory usage in bytes
 */
void get_analysis_performance(double *frames_per_second, double *realtime_factor,
                             size_t *memory_usage);

#endif /* AUDIO_ANALYSIS_H */

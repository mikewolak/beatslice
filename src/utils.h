/*
 * utils.h - Utility functions header
 * 
 * Mathematical, DSP, and general utility function declarations
 * for the beatslice audio analysis tool.
 */

#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>
#include <stdint.h>

/*----- Mathematical Utility Functions -----*/

/**
 * Check if a number is a power of two
 *
 * @param n Number to check
 * @return 1 if power of two, 0 otherwise
 */
int is_power_of_two(int n);

/**
 * Get current time with high precision
 *
 * @return Current time in seconds as double
 */
double get_current_time(void);

/**
 * Normalize a sample value
 *
 * @param sample Sample value to normalize
 * @param max_amplitude Maximum amplitude for normalization
 * @return Normalized sample value
 */
double normalize_sample(double sample, double max_amplitude);

/*----- DSP Utility Functions -----*/

/**
 * Hann window function
 *
 * @param n Sample index
 * @param N Total window size
 * @return Window coefficient value
 */
double hann_window(int n, int N);

/**
 * Generate Hamming window coefficients
 *
 * @param window Output array for window coefficients
 * @param length Window length
 */
void hamming_window(double *window, int length);

/**
 * Generate Blackman window coefficients
 *
 * @param window Output array for window coefficients  
 * @param length Window length
 */
void blackman_window(double *window, int length);

/**
 * Apply median filtering to a signal
 *
 * @param input Input signal array
 * @param output Output filtered signal array
 * @param length Length of signal arrays
 * @param filter_length Length of median filter
 */
void smooth_signal(const double *input, double *output, int length, int window_size);

/**
 * Calculate median value from buffer
 *
 * @param buffer Input buffer
 * @param length Buffer length
 * @return Median value
 */
double median_filter_sample(const double *buffer, int length);

/**
 * Interpolate peak location using quadratic interpolation
 *
 * @param left Left sample value
 * @param center Center sample value (peak)
 * @param right Right sample value
 * @return Sub-sample peak offset (-0.5 to 0.5)
 */
double interpolate_peak(double left, double center, double right);

/*----- Spectral Utility Functions -----*/

/**
 * Calculate spectral centroid
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
 * @param spectrum Magnitude spectrum array
 * @param length Length of spectrum array
 * @param sample_rate Audio sample rate
 * @param rolloff_point Rolloff percentage (e.g., 0.85)
 * @return Spectral rolloff frequency in Hz
 */
double calculate_spectral_rolloff(const double *spectrum, int length, 
                                 double sample_rate, double rolloff_point);

/**
 * Calculate spectral flatness measure
 *
 * @param spectrum Magnitude spectrum array
 * @param length Length of spectrum array
 * @return Spectral flatness value (0.0-1.0)
 */
double calculate_spectral_flatness(const double *spectrum, int length);

/**
 * Calculate RMS energy
 *
 * @param samples Audio sample array
 * @param length Number of samples
 * @return RMS energy value
 */
double calculate_rms(const float *samples, int length);

/**
 * Calculate zero crossing rate
 *
 * @param samples Audio sample array
 * @param length Number of samples
 * @return Zero crossing rate (0.0-1.0)
 */
double calculate_zero_crossing_rate(const float *samples, int length);

/**
 * Apply pre-emphasis filter
 *
 * @param samples Audio samples (modified in place)
 * @param length Number of samples
 * @param alpha Pre-emphasis coefficient (typically 0.95-0.97)
 */
void apply_pre_emphasis(float *samples, int length, double alpha);

/*----- Tempo and Beat Utility Functions -----*/

/**
 * Estimate tempo from beat positions
 *
 * @param beat_times Array of beat times in seconds
 * @param beat_count Number of beats
 * @param min_bpm Minimum BPM to consider
 * @param max_bpm Maximum BPM to consider
 * @return Estimated BPM value
 */
double estimate_tempo(const double *beat_times, int beat_count, 
                     double min_bpm, double max_bpm);

/*----- File and System Utilities -----*/

/**
 * Create output directory if it doesn't exist
 *
 * @param path Directory path to create
 * @return 0 on success, error code on failure
 */
int create_output_directory(const char *path);

/**
 * Get file size in bytes
 *
 * @param filename Path to file
 * @return File size in bytes, or -1 on error
 */
int64_t get_file_size(const char *filename);

/**
 * Check if filename has valid audio extension
 *
 * @param filename Filename to check
 * @return 1 if valid audio extension, 0 otherwise
 */
int has_valid_audio_extension(const char *filename);

/*----- Progress and Status Functions -----*/

/**
 * Progress callback for long operations
 *
 * @param percentage Progress percentage (0.0-100.0)
 * @param message Optional status message
 */
void progress_callback(double percentage, const char *message);

/*----- Error Handling Utilities -----*/

/* Note: beatslice_error_string is declared in beatslice.h */

/*----- Memory Management Utilities -----*/

/**
 * Safe memory allocation with error checking
 *
 * @param size Number of bytes to allocate
 * @return Pointer to allocated memory, or NULL on failure
 */
void* safe_malloc(size_t size);

/**
 * Safe memory reallocation with error checking
 *
 * @param ptr Existing pointer (can be NULL)
 * @param size New size in bytes
 * @return Pointer to reallocated memory, or NULL on failure
 */
void* safe_realloc(void *ptr, size_t size);

/**
 * Safe string duplication
 *
 * @param str String to duplicate
 * @return Pointer to duplicated string, or NULL on failure
 */
char* safe_strdup(const char *str);

/*----- Array and Buffer Utilities -----*/

/**
 * Find minimum value in array
 *
 * @param array Input array
 * @param length Array length
 * @return Minimum value
 */
double array_min(const double *array, int length);

/**
 * Find maximum value in array
 *
 * @param array Input array
 * @param length Array length
 * @return Maximum value
 */
double array_max(const double *array, int length);

/**
 * Calculate array mean
 *
 * @param array Input array
 * @param length Array length
 * @return Mean value
 */
double array_mean(const double *array, int length);

/**
 * Calculate array standard deviation
 *
 * @param array Input array
 * @param length Array length
 * @return Standard deviation
 */
double array_std(const double *array, int length);

/**
 * Zero-fill array
 *
 * @param array Array to zero
 * @param length Array length
 */
void array_zero(double *array, int length);

/**
 * Copy array contents
 *
 * @param dest Destination array
 * @param src Source array
 * @param length Number of elements to copy
 */
void array_copy(double *dest, const double *src, int length);

/*----- String Utilities -----*/

/**
 * Case-insensitive string comparison
 *
 * @param str1 First string
 * @param str2 Second string
 * @return 0 if equal, non-zero if different
 */
int string_compare_nocase(const char *str1, const char *str2);

/**
 * Extract file extension from filename
 *
 * @param filename Input filename
 * @return Pointer to extension (without dot), or NULL if none
 */
const char* get_file_extension(const char *filename);

/**
 * Extract base filename without extension
 *
 * @param filename Input filename
 * @param basename Output buffer for base name
 * @param buffer_size Size of output buffer
 * @return 0 on success, error code on failure
 */
int get_base_filename(const char *filename, char *basename, size_t buffer_size);

/*----- Debugging and Profiling -----*/

/**
 * Start timing measurement
 *
 * @return Timer handle for use with end_timing
 */
double start_timing(void);

/**
 * End timing measurement and get elapsed time
 *
 * @param start_time Timer handle from start_timing
 * @return Elapsed time in seconds
 */
double end_timing(double start_time);

/**
 * Print memory usage information
 */
void print_memory_usage(void);

/**
 * Log debug message (only in debug builds)
 *
 * @param format Printf-style format string
 * @param ... Format arguments
 */
void debug_log(const char *format, ...);

/*----- Mathematical Constants -----*/

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/*----- Utility Macros -----*/

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define CLAMP(x, min, max) MAX(min, MIN(max, x))
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

/* Safe division macro */
#define SAFE_DIVIDE(num, den, default) ((den) != 0.0 ? (num) / (den) : (default))

/* Decibel conversion macros */
#define LINEAR_TO_DB(x) (20.0 * log10(MAX(x, 1e-10)))
#define DB_TO_LINEAR(x) pow(10.0, (x) / 20.0)

#endif /* UTILS_H */

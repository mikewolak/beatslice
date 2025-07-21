/*
 * utils.c - Utility functions for beatslice
 * 
 * Mathematical, DSP, and general utility functions for audio analysis
 * and beat detection.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <errno.h>
#include "beatslice.h"
#include "utils.h"

/* Mathematical constants */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Static variables for progress reporting */
static double g_last_progress = -1.0;
static time_t g_start_time = 0;

/* Comparison function for qsort (used in median filter) */
static int compare_doubles(const void *a, const void *b) 
{
    double da = *(const double*)a;
    double db = *(const double*)b;
    
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

int is_power_of_two(int n) 
{
    return n > 0 && (n & (n - 1)) == 0;
}

double get_current_time(void) 
{
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return ts.tv_sec + ts.tv_nsec / 1e9;
    } else {
        /* Fallback to time() */
        return (double)time(NULL);
    }
}

void progress_callback(double percentage, const char *message) 
{
    /* Only update progress if it has changed significantly */
    if (fabs(percentage - g_last_progress) < 1.0) {
        return;
    }
    
    if (g_start_time == 0) {
        g_start_time = time(NULL);
    }
    
    time_t current_time = time(NULL);
    double elapsed = difftime(current_time, g_start_time);
    
    printf("\r%s: %.1f%%", message ? message : "Processing", percentage);
    
    if (percentage > 0 && elapsed > 1.0) {
        double eta = (elapsed / percentage) * (100.0 - percentage);
        if (eta < 60.0) {
            printf(" (ETA: %.0fs)", eta);
        } else {
            printf(" (ETA: %.0fm%.0fs)", floor(eta / 60.0), fmod(eta, 60.0));
        }
    }
    
    if (percentage >= 100.0) {
        printf(" - Done!\n");
        g_start_time = 0;
    } else {
        fflush(stdout);
    }
    
    g_last_progress = percentage;
}

int create_output_directory(const char *path) 
{
    if (!path) {
        return BEATSLICE_ERROR_INVALID_PARAMETER;
    }
    
    struct stat st;
    
    /* Check if directory already exists */
    if (stat(path, &st) == 0) {
        if (S_ISDIR(st.st_mode)) {
            return BEATSLICE_SUCCESS;
        } else {
            LOG_ERROR("Path exists but is not a directory: %s", path);
            return BEATSLICE_ERROR_INVALID_PARAMETER;
        }
    }
    
    /* Create directory */
    if (mkdir(path, 0755) != 0) {
        LOG_ERROR("Failed to create directory '%s': %s", path, strerror(errno));
        return BEATSLICE_ERROR_EXPORT_FAILED;
    }
    
    return BEATSLICE_SUCCESS;
}

double normalize_sample(double sample, double max_amplitude) 
{
    if (max_amplitude <= 0.0) {
        return 0.0;
    }
    
    return sample / max_amplitude;
}

double hann_window(int n, int N) 
{
    if (N <= 1) {
        return 1.0;
    }
    
    return 0.5 * (1.0 - cos(2.0 * M_PI * n / (N - 1)));
}

double median_filter_sample(const double *buffer, int length) 
{
    if (!buffer || length <= 0) {
        return 0.0;
    }
    
    if (length == 1) {
        return buffer[0];
    }
    
    /* Copy buffer for sorting */
    double *temp = malloc(sizeof(double) * length);
    if (!temp) {
        return buffer[length / 2];  /* Fallback to middle element */
    }
    
    memcpy(temp, buffer, sizeof(double) * length);
    
    /* Sort the buffer */
    qsort(temp, length, sizeof(double), compare_doubles);
    
    /* Return median */
    double result;
    if (length % 2 == 0) {
        /* Even length - average of two middle elements */
        result = (temp[length / 2 - 1] + temp[length / 2]) / 2.0;
    } else {
        /* Odd length - middle element */
        result = temp[length / 2];
    }
    
    free(temp);
    return result;
}

void smooth_signal(const double *input, double *output, int length, int window_size) 
{
    if (!input || !output || length <= 0 || window_size <= 0) {
        return;
    }
    
    int half_window = window_size / 2;
    
    for (int i = 0; i < length; i++) {
        double sum = 0.0;
        int count = 0;
        
        int start = fmax(0, i - half_window);
        int end = fmin(length - 1, i + half_window);
        
        for (int j = start; j <= end; j++) {
            sum += input[j];
            count++;
        }
        
        output[i] = (count > 0) ? sum / count : input[i];
    }
}

double calculate_spectral_centroid(const double *spectrum, int length, double sample_rate) 
{
    if (!spectrum || length <= 0 || sample_rate <= 0) {
        return 0.0;
    }
    
    double weighted_sum = 0.0;
    double magnitude_sum = 0.0;
    double freq_resolution = sample_rate / (2.0 * (length - 1));
    
    for (int i = 1; i < length; i++) {  /* Skip DC component */
        double frequency = i * freq_resolution;
        double magnitude = spectrum[i];
        
        weighted_sum += frequency * magnitude;
        magnitude_sum += magnitude;
    }
    
    return (magnitude_sum > 0.0) ? weighted_sum / magnitude_sum : 0.0;
}

double calculate_spectral_rolloff(const double *spectrum, int length, 
                                 double sample_rate, double rolloff_point) 
{
    if (!spectrum || length <= 0 || sample_rate <= 0 || 
        rolloff_point <= 0.0 || rolloff_point >= 1.0) {
        return 0.0;
    }
    
    /* Calculate total energy */
    double total_energy = 0.0;
    for (int i = 1; i < length; i++) {  /* Skip DC component */
        total_energy += spectrum[i];
    }
    
    if (total_energy <= 0.0) {
        return 0.0;
    }
    
    /* Find rolloff frequency */
    double target_energy = total_energy * rolloff_point;
    double cumulative_energy = 0.0;
    double freq_resolution = sample_rate / (2.0 * (length - 1));
    
    for (int i = 1; i < length; i++) {
        cumulative_energy += spectrum[i];
        
        if (cumulative_energy >= target_energy) {
            return i * freq_resolution;
        }
    }
    
    /* If not found, return Nyquist frequency */
    return sample_rate / 2.0;
}

double estimate_tempo(const double *beat_times, int beat_count, 
                     double min_bpm, double max_bpm) 
{
    if (!beat_times || beat_count < 2) {
        return 0.0;
    }
    
    /* Calculate inter-beat intervals */
    double *intervals = malloc(sizeof(double) * (beat_count - 1));
    if (!intervals) {
        return 0.0;
    }
    
    for (int i = 0; i < beat_count - 1; i++) {
        intervals[i] = beat_times[i + 1] - beat_times[i];
    }
    
    /* Find the most common interval using histogram approach */
    double min_interval = 60.0 / max_bpm;
    double max_interval = 60.0 / min_bpm;
    int bin_count = 100;
    double bin_width = (max_interval - min_interval) / bin_count;
    
    int *histogram = calloc(bin_count, sizeof(int));
    if (!histogram) {
        free(intervals);
        return 0.0;
    }
    
    /* Populate histogram */
    for (int i = 0; i < beat_count - 1; i++) {
        if (intervals[i] >= min_interval && intervals[i] <= max_interval) {
            int bin = (int)((intervals[i] - min_interval) / bin_width);
            if (bin >= 0 && bin < bin_count) {
                histogram[bin]++;
            }
        }
    }
    
    /* Find the bin with maximum count */
    int max_bin = 0;
    int max_count = histogram[0];
    
    for (int i = 1; i < bin_count; i++) {
        if (histogram[i] > max_count) {
            max_count = histogram[i];
            max_bin = i;
        }
    }
    
    /* Calculate BPM from the most common interval */
    double estimated_interval = min_interval + (max_bin + 0.5) * bin_width;
    double estimated_bpm = 60.0 / estimated_interval;
    
    /* Validate result */
    if (estimated_bpm < min_bpm || estimated_bpm > max_bpm) {
        /* Fallback to simple average */
        double total_time = beat_times[beat_count - 1] - beat_times[0];
        estimated_bpm = (beat_count - 1) * 60.0 / total_time;
    }
    
    free(intervals);
    free(histogram);
    
    return estimated_bpm;
}

const char* beatslice_error_string(BeatsliceError error) 
{
    switch (error) {
        case BEATSLICE_SUCCESS:
            return "Success";
        case BEATSLICE_ERROR_FILE_NOT_FOUND:
            return "File not found or cannot be opened";
        case BEATSLICE_ERROR_UNSUPPORTED_FORMAT:
            return "Unsupported audio format";
        case BEATSLICE_ERROR_MEMORY_ALLOCATION:
            return "Memory allocation failed";
        case BEATSLICE_ERROR_INVALID_PARAMETER:
            return "Invalid parameter";
        case BEATSLICE_ERROR_ANALYSIS_FAILED:
            return "Audio analysis failed";
        case BEATSLICE_ERROR_EXPORT_FAILED:
            return "Export operation failed";
        default:
            return "Unknown error";
    }
}

/* Advanced mathematical functions for enhanced analysis */

double calculate_rms(const float *samples, int length) 
{
    if (!samples || length <= 0) {
        return 0.0;
    }
    
    double sum_squares = 0.0;
    for (int i = 0; i < length; i++) {
        double sample = samples[i];
        sum_squares += sample * sample;
    }
    
    return sqrt(sum_squares / length);
}

double calculate_zero_crossing_rate(const float *samples, int length) 
{
    if (!samples || length <= 1) {
        return 0.0;
    }
    
    int crossings = 0;
    for (int i = 1; i < length; i++) {
        if ((samples[i-1] >= 0.0) != (samples[i] >= 0.0)) {
            crossings++;
        }
    }
    
    return (double)crossings / (length - 1);
}

void apply_pre_emphasis(float *samples, int length, double alpha) 
{
    if (!samples || length <= 1 || alpha <= 0.0 || alpha >= 1.0) {
        return;
    }
    
    for (int i = length - 1; i > 0; i--) {
        samples[i] = samples[i] - alpha * samples[i - 1];
    }
}

double calculate_spectral_flatness(const double *spectrum, int length) 
{
    if (!spectrum || length <= 0) {
        return 0.0;
    }
    
    double geometric_mean = 0.0;
    double arithmetic_mean = 0.0;
    int valid_bins = 0;
    
    for (int i = 1; i < length; i++) {  /* Skip DC component */
        if (spectrum[i] > 0.0) {
            geometric_mean += log(spectrum[i]);
            arithmetic_mean += spectrum[i];
            valid_bins++;
        }
    }
    
    if (valid_bins == 0) {
        return 0.0;
    }
    
    geometric_mean = exp(geometric_mean / valid_bins);
    arithmetic_mean /= valid_bins;
    
    return (arithmetic_mean > 0.0) ? geometric_mean / arithmetic_mean : 0.0;
}

void hamming_window(double *window, int length) 
{
    if (!window || length <= 0) {
        return;
    }
    
    for (int i = 0; i < length; i++) {
        window[i] = 0.54 - 0.46 * cos(2.0 * M_PI * i / (length - 1));
    }
}

void blackman_window(double *window, int length) 
{
    if (!window || length <= 0) {
        return;
    }
    
    for (int i = 0; i < length; i++) {
        double factor = 2.0 * M_PI * i / (length - 1);
        window[i] = 0.42 - 0.5 * cos(factor) + 0.08 * cos(2.0 * factor);
    }
}

double interpolate_peak(double left, double center, double right) 
{
    /* Quadratic interpolation for sub-sample peak estimation */
    double denominator = 2.0 * (2.0 * center - left - right);
    
    if (fabs(denominator) < 1e-10) {
        return 0.0;  /* No interpolation possible */
    }
    
    return (right - left) / denominator;
}

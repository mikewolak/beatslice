/* PART 1 START */
/* src/audio_analysis.c - Audio analysis implementation with MULTI-BAND onset detection */
#include "beatslice.h"
#include "audio_analysis.h"
#include "utils.h"
#include <math.h>
#include <string.h>

/* Define M_PI if not already defined (ANSI C90 doesn't guarantee it) */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Internal spectral frame buffer - adapted to new structures */
static SpectralFrame g_frames[MAX_ANALYSIS_FRAMES];
static int g_frame_count = 0;

/* Multi-band onset detection arrays */
static double g_low_envelope[MAX_ANALYSIS_FRAMES];
static double g_mid_envelope[MAX_ANALYSIS_FRAMES];
static double g_high_envelope[MAX_ANALYSIS_FRAMES];

int detect_beats_multiband_optimized(const AudioData *audio, const BeatsliceConfig *config,
                                    AnalysisResults *results);

/* Global analysis state - adapted for FFTW3 */
static struct {
    fftw_plan fft_plan;
    double *fft_input;
    fftw_complex *fft_output;
    double *spectrum;
    double *energy_envelope;
    int fft_initialized;
    int fft_size;
    int hop_size;
    int sample_rate;
} g_analysis_state = {0};

/**
 * Initialize the FFT for spectral analysis - ADAPTED FOR FFTW3
 */
int analysis_init(const BeatsliceConfig *config)
{
    LOG_DEBUG("analysis_init: Starting with FFT size %d", config->fft_size);
    
    /* Check if already initialized */
    if (g_analysis_state.fft_initialized) {
        LOG_DEBUG("analysis_init: Already initialized");
        analysis_cleanup();
    }
    
    g_analysis_state.fft_size = config->fft_size;
    g_analysis_state.hop_size = config->hop_size;
    
    /* Allocate FFT buffers - FFTW3 style */
    g_analysis_state.fft_input = fftw_malloc(sizeof(double) * config->fft_size);
    g_analysis_state.fft_output = fftw_malloc(sizeof(fftw_complex) * (config->fft_size/2 + 1));
    
    if (!g_analysis_state.fft_input || !g_analysis_state.fft_output) {
        LOG_ERROR("analysis_init: Failed to allocate FFT buffers");
        return BEATSLICE_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Create FFTW3 plan - ADAPTED */
    g_analysis_state.fft_plan = fftw_plan_dft_r2c_1d(config->fft_size, 
                                                     g_analysis_state.fft_input,
                                                     g_analysis_state.fft_output, 
                                                     FFTW_ESTIMATE);
    if (!g_analysis_state.fft_plan) {
        LOG_ERROR("analysis_init: Failed to create FFTW plan");
        fftw_free(g_analysis_state.fft_input);
        fftw_free(g_analysis_state.fft_output);
        return BEATSLICE_ERROR_ANALYSIS_FAILED;
    }
    
    /* Allocate spectrum and energy envelope buffers */
    g_analysis_state.spectrum = malloc(sizeof(double) * (config->fft_size/2));
    g_analysis_state.energy_envelope = malloc(sizeof(double) * MAX_ANALYSIS_FRAMES);
    
    if (!g_analysis_state.spectrum || !g_analysis_state.energy_envelope) {
        LOG_ERROR("analysis_init: Failed to allocate analysis buffers");
        fftw_destroy_plan(g_analysis_state.fft_plan);
        fftw_free(g_analysis_state.fft_input);
        fftw_free(g_analysis_state.fft_output);
        return BEATSLICE_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Mark as initialized */
    g_analysis_state.fft_initialized = 1;
    LOG_DEBUG("analysis_init: Successfully initialized FFT with size %d", config->fft_size);
    
    return BEATSLICE_SUCCESS;
}

/**
 * Clean up FFT resources - ADAPTED FOR FFTW3
 */
void analysis_cleanup(void)
{
    LOG_DEBUG("analysis_cleanup: Starting");
    
    if (!g_analysis_state.fft_initialized) {
        LOG_DEBUG("analysis_cleanup: FFT not initialized");
        return;
    }
    
    /* Free FFTW3 resources */
    if (g_analysis_state.fft_plan) {
        fftw_destroy_plan(g_analysis_state.fft_plan);
        g_analysis_state.fft_plan = NULL;
    }
    
    if (g_analysis_state.fft_input) {
        fftw_free(g_analysis_state.fft_input);
        g_analysis_state.fft_input = NULL;
    }
    
    if (g_analysis_state.fft_output) {
        fftw_free(g_analysis_state.fft_output);
        g_analysis_state.fft_output = NULL;
    }
    
    free(g_analysis_state.spectrum);
    free(g_analysis_state.energy_envelope);
    g_analysis_state.spectrum = NULL;
    g_analysis_state.energy_envelope = NULL;
    
    g_analysis_state.fft_initialized = 0;
    LOG_DEBUG("analysis_cleanup: Successfully cleaned up FFT resources");
}

/**
 * Apply a Hann window to an audio frame - UNCHANGED
 */
static void apply_hann_window(double *buffer, int size)
{
    int i;
    double multiplier;
    
    for (i = 0; i < size; i++) {
        /* Hann window: 0.5 * (1 - cos(2*pi*n/(N-1))) */
        multiplier = 0.5 * (1.0 - cos(2.0 * M_PI * i / (size - 1)));
        buffer[i] *= multiplier;
    }
}
/* PART 1 END */

/* PART 2 START */
/**
 * Extract samples from audio data - PROPER FLOAT HANDLING
 * 
 * Keep libsndfile normalized float data as-is for modern, clean processing
 */
static void extract_samples(const AudioData *audio, int offset, double *output, int size)
{
    int i;
    
    /* Zero output buffer */
    memset(output, 0, size * sizeof(double));
    
    /* Check offset bounds */
    if (offset >= audio->frame_count) {
        return;
    }
    
    /* Limit size if it would exceed sample count */
    if (offset + size > audio->frame_count) {
        size = audio->frame_count - offset;
    }
    
    /* Extract samples - KEEP FLOAT RANGE (-1.0 to +1.0) */
    if (audio->channels == 2) {
        /* Stereo - average the channels */
        for (i = 0; i < size; i++) {
            if (offset + i < audio->frame_count) {
                /* Average left and right channels - keep float range */
                double left = audio->samples[(offset + i) * 2];
                double right = audio->samples[(offset + i) * 2 + 1];
                output[i] = (left + right) / 2.0;
            }
        }
    } else {
        /* Mono - direct copy */
        for (i = 0; i < size; i++) {
            if (offset + i < audio->frame_count) {
                output[i] = audio->samples[offset + i];
            }
        }
    }
}

/**
 * Calculate RMS energy for audio data - PROPER FLOAT HANDLING
 */
static double calculate_frame_rms(const AudioData *audio, int offset, int window_size)
{
    int i;
    double sum = 0.0;
    double sample;
    
    for (i = 0; i < window_size; i++) {
        if (offset + i < audio->frame_count) {
            if (audio->channels == 2) {
                /* Average stereo channels */
                double left = audio->samples[(offset + i) * 2];
                double right = audio->samples[(offset + i) * 2 + 1];
                sample = (left + right) / 2.0;
            } else {
                sample = audio->samples[offset + i];
            }
            sum += sample * sample;
        }
    }
    
    return sqrt(sum / window_size);
}

/**
 * Count zero crossings in audio data - PROPER FLOAT HANDLING
 */
static int count_zero_crossings(const AudioData *audio, int offset, int window_size)
{
    int i;
    int crossings = 0;
    int prev_sign, curr_sign;
    double sample, prev_sample;
    
    if (window_size < 2) return 0;
    
    /* Get first sample */
    if (audio->channels == 2) {
        double left = audio->samples[offset * 2];
        double right = audio->samples[offset * 2 + 1];
        prev_sample = (left + right) / 2.0;
    } else {
        prev_sample = audio->samples[offset];
    }
    prev_sign = (prev_sample > 0.0) ? 1 : ((prev_sample < 0.0) ? -1 : 0);
    
    for (i = 1; i < window_size; i++) {
        if (offset + i < audio->frame_count) {
            if (audio->channels == 2) {
                double left = audio->samples[(offset + i) * 2];
                double right = audio->samples[(offset + i) * 2 + 1];
                sample = (left + right) / 2.0;
            } else {
                sample = audio->samples[offset + i];
            }
            
            curr_sign = (sample > 0.0) ? 1 : ((sample < 0.0) ? -1 : 0);
            
            /* Count only when sign changes from positive to negative or vice versa */
            if (prev_sign != 0 && curr_sign != 0 && prev_sign != curr_sign) {
                crossings++;
            }
            
            prev_sign = curr_sign;
        }
    }
    
    return crossings;
}
/* PART 2 END */

/* PART 3 START */
/**
 * Extract energy in specific frequency band - ADAPTED FOR FFTW3
 */
static double get_band_energy(fftw_complex *fft_out, int fft_size, int sample_rate, 
                              double low_freq, double high_freq)
{
    int low_bin, high_bin;
    int i;
    double energy = 0.0;
    double real, imag, magnitude;
    double bin_width = (double)sample_rate / fft_size;
    
    /* Calculate bin indices for frequency range */
    low_bin = (int)(low_freq / bin_width);
    high_bin = (int)(high_freq / bin_width);
    
    /* Clamp to valid range */
    if (low_bin < 1) low_bin = 1;  /* Skip DC component (bin 0) */
    if (high_bin >= fft_size/2) high_bin = fft_size/2 - 1;
    
    /* Sum energy in the bins - FFTW3 format */
    for (i = low_bin; i <= high_bin; i++) {
        /* FFTW3 complex format: fft_out[i][0] = real, fft_out[i][1] = imag */
        real = fft_out[i][0];
        imag = fft_out[i][1];
        magnitude = sqrt(real*real + imag*imag);
        energy += magnitude;
    }
    
    /* Normalize by the number of bins */
    if (high_bin >= low_bin) {
        energy /= (high_bin - low_bin + 1);
    }
    
    return energy;
}

/**
 * Perform spectral analysis on loaded audio - MULTI-BAND ONSET DETECTION
 */
int perform_beat_analysis(const AudioData *audio, const BeatsliceConfig *config, 
                         AnalysisResults *results)
{
    int i;
    int hop_samples;
    int window_offset;
    double max_low_energy = 0.0, max_mid_energy = 0.0, max_high_energy = 0.0;
    
    LOG_DEBUG("perform_beat_analysis: Starting MULTI-BAND analysis with FFT size %d", config->fft_size);
    
    /* Check if we have audio and FFT is initialized */
    if (!audio->samples || !g_analysis_state.fft_initialized) {
        LOG_ERROR("perform_beat_analysis: No audio or FFT not initialized");
        return BEATSLICE_ERROR_ANALYSIS_FAILED;
    }
    
    /* Set sample rate for global state */
    g_analysis_state.sample_rate = audio->sample_rate;
    
    /* Calculate hop size in samples */
    hop_samples = config->hop_size;
    
    LOG_DEBUG("perform_beat_analysis: Starting analysis with hop size %d", hop_samples);
    LOG_DEBUG("Frequency bands: LOW(40-250Hz), MID(250-2000Hz), HIGH(2000-%.0fHz)", 
              audio->sample_rate / 2.0);
    
    /* Calculate number of frames */
    int spectrum_frames = (audio->frame_count - config->fft_size) / hop_samples + 1;
    
    /* Limit to maximum frames */
    if (spectrum_frames > MAX_ANALYSIS_FRAMES) {
        spectrum_frames = MAX_ANALYSIS_FRAMES;
    }
    
    LOG_DEBUG("perform_beat_analysis: Processing %d frames", spectrum_frames);
    
    /* Process each frame */
    for (i = 0; i < spectrum_frames; i++) {
        /* Calculate window offset */
        window_offset = i * hop_samples;
        
        /* Extract samples - KEEP FLOAT RANGE */
        extract_samples(audio, window_offset, g_analysis_state.fft_input, config->fft_size);
        
        /* Apply window function */
        apply_hann_window(g_analysis_state.fft_input, config->fft_size);
        
        /* Perform FFT - FFTW3 style */
        fftw_execute(g_analysis_state.fft_plan);
        
        /* Calculate time position */
        g_frames[i].time_position = (double)window_offset / audio->sample_rate;
        
        /* Calculate RMS for the frame */
        g_frames[i].rms_energy = calculate_frame_rms(audio, window_offset, config->fft_size);
        
        /* Count zero crossings if enabled */
        if (config->use_zero_crossings) {
            g_frames[i].zero_crossing_rate = (double)count_zero_crossings(audio, window_offset, config->fft_size) / config->fft_size;
        } else {
            g_frames[i].zero_crossing_rate = 0.0;
        }
        
        /* Calculate band energies - ALL THREE BANDS */
        g_frames[i].low_energy = get_band_energy(g_analysis_state.fft_output,
                                                config->fft_size,
                                                audio->sample_rate,
                                                40.0, 250.0);  /* Bass/kick range */
        
        g_frames[i].mid_energy = get_band_energy(g_analysis_state.fft_output,
                                               config->fft_size,
                                               audio->sample_rate,
                                               250.0, 2000.0);  /* Snare/vocal range */
        
        g_frames[i].high_energy = get_band_energy(g_analysis_state.fft_output,
                                                config->fft_size,
                                                audio->sample_rate,
                                                2000.0, audio->sample_rate / 2.0);  /* Hi-hat/cymbal range */
        
        /* Store in separate band envelopes */
        g_low_envelope[i] = g_frames[i].low_energy;
        g_mid_envelope[i] = g_frames[i].mid_energy;
        g_high_envelope[i] = g_frames[i].high_energy;
        
        /* Track maximum energies for normalization */
        if (g_low_envelope[i] > max_low_energy) max_low_energy = g_low_envelope[i];
        if (g_mid_envelope[i] > max_mid_energy) max_mid_energy = g_mid_envelope[i];
        if (g_high_envelope[i] > max_high_energy) max_high_energy = g_high_envelope[i];
        
        /* Update progress */
        if (config->verbose && i % (spectrum_frames / 10) == 0) {
            progress_callback((double)(i * 100) / spectrum_frames, "Computing multi-band spectral features");
        }
    }
    
    /* Store frame count */
    g_frame_count = spectrum_frames;
    
    LOG_DEBUG("Multi-band energy maxima: LOW=%.6f, MID=%.6f, HIGH=%.6f", 
              max_low_energy, max_mid_energy, max_high_energy);
    
    /* Normalize each band envelope separately */
    if (max_low_energy > 0.0) {
        for (i = 0; i < spectrum_frames; i++) {
            g_low_envelope[i] /= max_low_energy;
        }
        LOG_DEBUG("Normalized LOW frequency envelope");
    }
    
    if (max_mid_energy > 0.0) {
        for (i = 0; i < spectrum_frames; i++) {
            g_mid_envelope[i] /= max_mid_energy;
        }
        LOG_DEBUG("Normalized MID frequency envelope");
    }
    
    if (max_high_energy > 0.0) {
        for (i = 0; i < spectrum_frames; i++) {
            g_high_envelope[i] /= max_high_energy;
        }
        LOG_DEBUG("Normalized HIGH frequency envelope");
    }
    
    /* Check if we have any energy in any band */
    if (max_low_energy <= 0.0 && max_mid_energy <= 0.0 && max_high_energy <= 0.0) {
        LOG_ERROR("perform_beat_analysis: No energy detected in any frequency band");
        return BEATSLICE_ERROR_ANALYSIS_FAILED;
    }
    
    /* Now perform multi-band beat detection */
    return detect_beats_multiband_optimized(audio, config, results);
}
/* PART 3 END */

/* PART 4 START */
/**
 * Create adaptive threshold for multi-band peak detection
 */
static void create_adaptive_threshold_multiband(double *low_env, double *mid_env, double *high_env,
                                               int size, int window_size, double threshold_factor,
                                               double *low_thresh, double *mid_thresh, double *high_thresh)
{
    int i, j;
    int half_window = window_size / 2;
    
    for (i = 0; i < size; i++) {
        double low_sum = 0.0, mid_sum = 0.0, high_sum = 0.0;
        int count = 0;
        
        /* Calculate average energy in the window for each band */
        for (j = i - half_window; j <= i + half_window; j++) {
            if (j >= 0 && j < size) {
                low_sum += low_env[j];
                mid_sum += mid_env[j];
                high_sum += high_env[j];
                count++;
            }
        }
        
        /* Calculate thresholds - different factors for different bands */
        if (count > 0) {
            low_thresh[i] = (low_sum / count) * threshold_factor;
            mid_thresh[i] = (mid_sum / count) * (threshold_factor * 0.8);  /* Slightly lower for mid */
            high_thresh[i] = (high_sum / count) * (threshold_factor * 0.6); /* Lower for high freq */
        } else {
            low_thresh[i] = mid_thresh[i] = high_thresh[i] = 0.0;
        }
    }
}

/**
 * Multi-band beat detection with detailed debug output
 */
int detect_beats_multiband_optimized(const AudioData *audio, const BeatsliceConfig *config,
                                    AnalysisResults *results)
{
    int i;
    int peak_indices[MAX_BEATS];
    char peak_bands[MAX_BEATS];  /* 'L', 'M', or 'H' for each peak */
    int peak_count;
    double low_threshold, mid_threshold, high_threshold;
    double *adaptive_low_thresh = NULL, *adaptive_mid_thresh = NULL, *adaptive_high_thresh = NULL;
    int min_peak_distance;
    double peak_threshold_factor;
    double min_peak_distance_sec = 0.3;  /* 300ms minimum like original */
    
    LOG_DEBUG("detect_beats_multiband_optimized: Starting MULTI-BAND detection");
    
    /* Check if we have spectral frames */
    if (g_frame_count == 0) {
        LOG_ERROR("detect_beats_multiband_optimized: No spectral frames available");
        return BEATSLICE_ERROR_ANALYSIS_FAILED;
    }
    
    /* Configure beat detector - MULTI-BAND OPTIMIZED THRESHOLDS */
    peak_threshold_factor = 0.3 + (config->sensitivity * 0.5);  /* Range: 0.3 to 0.8 */
    
    LOG_DEBUG("Multi-band detection: Threshold factor = %.3f", peak_threshold_factor);
    
    /* Calculate minimum distance between peaks in frames */
    min_peak_distance = (int)(min_peak_distance_sec * audio->sample_rate / config->hop_size);
    if (min_peak_distance < 1) min_peak_distance = 1;
    
    LOG_DEBUG("Multi-band detection: Minimum peak distance = %d frames", min_peak_distance);
    
    /* Create adaptive thresholds if enabled */
    if (config->use_adaptive_threshold) {
        adaptive_low_thresh = malloc(sizeof(double) * g_frame_count);
        adaptive_mid_thresh = malloc(sizeof(double) * g_frame_count);
        adaptive_high_thresh = malloc(sizeof(double) * g_frame_count);
        
        if (adaptive_low_thresh && adaptive_mid_thresh && adaptive_high_thresh) {
            create_adaptive_threshold_multiband(g_low_envelope, g_mid_envelope, g_high_envelope,
                                               g_frame_count, min_peak_distance * 2, 
                                               peak_threshold_factor,
                                               adaptive_low_thresh, adaptive_mid_thresh, adaptive_high_thresh);
            
            LOG_DEBUG("Multi-band adaptive thresholds created");
            
            /* Find peaks using adaptive thresholds */
            peak_count = 0;
            memset(peak_bands, 0, sizeof(peak_bands));
            
            for (i = 1; i < g_frame_count - 1; i++) {
                double low_val = g_low_envelope[i];
                double mid_val = g_mid_envelope[i];
                double high_val = g_high_envelope[i];
                
                char triggered_band = 0;
                double max_val = 0.0;
                
                /* Check each band against its adaptive threshold */
                if (low_val > adaptive_low_thresh[i] && 
                    low_val > g_low_envelope[i-1] && low_val >= g_low_envelope[i+1]) {
                    if (low_val > max_val) {
                        max_val = low_val;
                        triggered_band = 'L';
                    }
                }
                
                if (mid_val > adaptive_mid_thresh[i] && 
                    mid_val > g_mid_envelope[i-1] && mid_val >= g_mid_envelope[i+1]) {
                    if (mid_val > max_val) {
                        max_val = mid_val;
                        triggered_band = 'M';
                    }
                }
                
                if (high_val > adaptive_high_thresh[i] && 
                    high_val > g_high_envelope[i-1] && high_val >= g_high_envelope[i+1]) {
                    if (high_val > max_val) {
                        max_val = high_val;
                        triggered_band = 'H';
                    }
                }
                
                if (triggered_band && peak_count < MAX_BEATS) {
                    /* Check minimum distance constraint */
                    int valid = 1;
                    for (int j = 0; j < peak_count; j++) {
                        if (abs(i - peak_indices[j]) < min_peak_distance) {
                            valid = 0;
                            break;
                        }
                    }
                    
                    if (valid) {
                        peak_indices[peak_count] = i;
                        peak_bands[peak_count] = triggered_band;
                        
                        /* DEBUG OUTPUT - Show which band triggered */
                        double time_pos = (double)(i * config->hop_size) / audio->sample_rate;
                        const char *band_name = (triggered_band == 'L') ? "LOW" : 
                                               (triggered_band == 'M') ? "MID" : "HIGH";
                        
                        if (config->verbose) {
                            LOG_INFO("ONSET DETECTED: %.3fs in %s frequency band (%.3f energy)", 
                                    time_pos, band_name, max_val);
                        }
                        
                        peak_count++;
                    }
                }
            }
            
            free(adaptive_low_thresh);
            free(adaptive_mid_thresh);
            free(adaptive_high_thresh);
        } else {
            /* Fallback to fixed thresholds */
            LOG_WARNING("Failed to allocate adaptive threshold arrays, using fixed thresholds");
            goto fixed_threshold_fallback;
        }
    } else {
fixed_threshold_fallback:
        /* Use fixed thresholds - MULTI-BAND OPTIMIZED */
        low_threshold = 0.15 + (config->sensitivity * 0.25);   /* Range: 0.15 to 0.4 */
        mid_threshold = 0.12 + (config->sensitivity * 0.23);   /* Slightly lower for mid */
        high_threshold = 0.10 + (config->sensitivity * 0.20);  /* Lower for high freq */
        
        LOG_DEBUG("Fixed thresholds: LOW=%.3f, MID=%.3f, HIGH=%.3f", 
                  low_threshold, mid_threshold, high_threshold);
        
        /* Find peaks in each band separately */
        peak_count = 0;
        memset(peak_bands, 0, sizeof(peak_bands));
        
        for (i = 1; i < g_frame_count - 1; i++) {
            char triggered_band = 0;
            double max_val = 0.0;
            
            /* Check LOW band */
            if (g_low_envelope[i] > low_threshold && 
                g_low_envelope[i] > g_low_envelope[i-1] && 
                g_low_envelope[i] >= g_low_envelope[i+1]) {
                max_val = g_low_envelope[i];
                triggered_band = 'L';
            }
            
            /* Check MID band */
            if (g_mid_envelope[i] > mid_threshold && 
                g_mid_envelope[i] > g_mid_envelope[i-1] && 
                g_mid_envelope[i] >= g_mid_envelope[i+1] &&
                g_mid_envelope[i] > max_val) {
                max_val = g_mid_envelope[i];
                triggered_band = 'M';
            }
            
            /* Check HIGH band */
            if (g_high_envelope[i] > high_threshold && 
                g_high_envelope[i] > g_high_envelope[i-1] && 
                g_high_envelope[i] >= g_high_envelope[i+1] &&
                g_high_envelope[i] > max_val) {
                max_val = g_high_envelope[i];
                triggered_band = 'H';
            }
            
            if (triggered_band && peak_count < MAX_BEATS) {
                /* Check minimum distance constraint */
                int valid = 1;
                for (int j = 0; j < peak_count; j++) {
                    if (abs(i - peak_indices[j]) < min_peak_distance) {
                        valid = 0;
                        break;
                    }
                }
                
                if (valid) {
                    peak_indices[peak_count] = i;
                    peak_bands[peak_count] = triggered_band;
                    
                    /* DEBUG OUTPUT - Show which band triggered */
                    double time_pos = (double)(i * config->hop_size) / audio->sample_rate;
                    const char *band_name = (triggered_band == 'L') ? "LOW" : 
                                           (triggered_band == 'M') ? "MID" : "HIGH";
                    
                    if (config->verbose) {
                        LOG_INFO("ONSET DETECTED: %.3fs in %s frequency band (%.3f energy)", 
                                time_pos, band_name, max_val);
                    }
                    
                    peak_count++;
                }
            }
        }
    }
    
    LOG_DEBUG("Multi-band detection found %d peaks", peak_count);
    
    /* Count peaks by band for summary */
    int low_count = 0, mid_count = 0, high_count = 0;
    for (i = 0; i < peak_count; i++) {
        switch (peak_bands[i]) {
            case 'L': low_count++; break;
            case 'M': mid_count++; break;
            case 'H': high_count++; break;
        }
    }
    
    LOG_INFO("Multi-band summary: %d LOW, %d MID, %d HIGH frequency onsets detected", 
             low_count, mid_count, high_count);
    
    /* Allocate results arrays */
    results->beat_times = malloc(sizeof(double) * peak_count);
    results->confidence = malloc(sizeof(double) * peak_count);
    
    if (!results->beat_times || !results->confidence) {
        return BEATSLICE_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Convert peak indices to beat positions with band information */
    double hop_time = (double)config->hop_size / audio->sample_rate;
    for (i = 0; i < peak_count; i++) {
        results->beat_times[i] = g_frames[peak_indices[i]].time_position;
        
        /* Set confidence based on the band that triggered */
        switch (peak_bands[i]) {
            case 'L': 
                results->confidence[i] = g_low_envelope[peak_indices[i]]; 
                break;
            case 'M': 
                results->confidence[i] = g_mid_envelope[peak_indices[i]]; 
                break;
            case 'H': 
                results->confidence[i] = g_high_envelope[peak_indices[i]]; 
                break;
            default:
                results->confidence[i] = 0.5;  /* Fallback */
        }
    }
    
    results->beat_count = peak_count;
    
    /* Calculate BPM - SAME AS ORIGINAL */
    if (peak_count > 1) {
        double total_time = results->beat_times[peak_count-1] - results->beat_times[0];
        results->estimated_bpm = (peak_count - 1) * 60.0 / total_time;
    } else {
        results->estimated_bpm = 0.0;
    }
    
    results->tempo_confidence = (peak_count > 2) ? 0.8 : 0.0;
    
    /* Store analysis data for export - use combined envelope */
    results->energy_envelope = malloc(sizeof(double) * g_frame_count);
    if (results->energy_envelope) {
        /* Create combined envelope for export (maximum of all bands) */
        for (i = 0; i < g_frame_count; i++) {
            double max_env = g_low_envelope[i];
            if (g_mid_envelope[i] > max_env) max_env = g_mid_envelope[i];
            if (g_high_envelope[i] > max_env) max_env = g_high_envelope[i];
            results->energy_envelope[i] = max_env;
        }
    }
    
    results->frame_count = g_frame_count;
    results->analysis_hop_time = hop_time;
    
    LOG_DEBUG("Multi-band beat detection complete: %d beats, estimated BPM: %.1f", 
              results->beat_count, results->estimated_bpm);
    
    return BEATSLICE_SUCCESS;
}
/* PART 4 END */

/* Stub implementations for compatibility */
void configure_sensitivity(double sensitivity, const BeatsliceConfig *config,
                          EnhancedBeatDetectionConfig *beat_config) 
{
    (void)sensitivity; (void)config; (void)beat_config; /* unused */
}

int compute_spectral_features(const AudioData *audio, const BeatsliceConfig *config,
                             SpectralFrame *frames, int *frame_count) 
{
    (void)audio; (void)config; (void)frames; (void)frame_count; /* unused */
    return BEATSLICE_SUCCESS;
}

int detect_onsets(const SpectralFrame *frames, int frame_count,
                 const EnhancedBeatDetectionConfig *config, double *onset_function) 
{
    (void)frames; (void)frame_count; (void)config; (void)onset_function; /* unused */
    return BEATSLICE_SUCCESS;
}

int detect_beats_enhanced(const double *onset_function, int function_length,
                         double hop_time, const EnhancedBeatDetectionConfig *config,
                         AnalysisResults *results) 
{
    (void)onset_function; (void)function_length; (void)hop_time; (void)config; (void)results;
    return BEATSLICE_SUCCESS;
}

double calculate_onset_strength(const SpectralFrame *current, 
                               const SpectralFrame *previous) 
{
    (void)current; (void)previous;
    return 0.0;
}

double calculate_spectral_flux(const double *current_spectrum, 
                              const double *previous_spectrum, int length) 
{
    (void)current_spectrum; (void)previous_spectrum; (void)length;
    return 0.0;
}

double calculate_high_frequency_content(const double *spectrum, int length) 
{
    (void)spectrum; (void)length;
    return 0.0;
}

void adaptive_threshold(const double *signal, int length, double factor,
                       double *threshold) 
{
    (void)signal; (void)length; (void)factor; (void)threshold;
}

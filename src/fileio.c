/*
 * fileio.c - Audio file I/O operations using libsndfile
 * 
 * Handles loading audio files in various formats and exporting beat slices
 * with support for multiple output formats and ATTACK-OPTIMIZED slicing.
 * 
 * ENHANCED: Attack-focused boundary alignment prevents audio clicks/pops
 * while preserving drum machine timing and immediate attack response.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include <libgen.h>
#include <math.h>
#include "beatslice.h"
#include "fileio.h"

/* Attack detection and zero crossing parameters */
#define MAX_ZERO_CROSSING_SEARCH_MS  5.0    /* Reduced to 5ms for tight alignment */
#define MIN_ZERO_CROSSING_AMPLITUDE  0.001  /* Minimum amplitude to consider as "zero" */
#define MAX_PRE_ATTACK_MS            3.0    /* Maximum 3ms of pre-attack data allowed */
#define ATTACK_DETECTION_WINDOW_MS   10.0   /* 10ms window to scan for attack */
#define ATTACK_ENERGY_THRESHOLD      2.0    /* Energy must increase by this factor for attack */

/* Output format mappings */
typedef struct {
    const char *extension;
    int sf_format;
    int sf_subtype;
} OutputFormat;

static const OutputFormat output_formats[] = {
    {"wav",  SF_FORMAT_WAV,  SF_FORMAT_PCM_16},
    {"aiff", SF_FORMAT_AIFF, SF_FORMAT_PCM_16},
    {"flac", SF_FORMAT_FLAC, SF_FORMAT_PCM_16},
    {NULL, 0, 0}
};

/* Private function prototypes */
static int get_output_format(const char *format_name, int *sf_format, int *sf_subtype);
static int ensure_output_directory(const char *filepath);
static void normalize_audio_slice(float *samples, sf_count_t frame_count, 
                                 int channels);
static sf_count_t find_nearest_zero_crossing(const AudioData *audio, sf_count_t target_frame,
                                            sf_count_t max_search_frames, int prefer_direction);
static sf_count_t find_attack_point(const AudioData *audio, sf_count_t target_frame, 
                                   sf_count_t search_window_frames);
static sf_count_t find_optimal_slice_boundary(const AudioData *audio, sf_count_t beat_frame, 
                                             int is_start_boundary);

/**
 * Find the nearest zero crossing to a target frame position
 *
 * @param audio Audio data to search in
 * @param target_frame Target frame position
 * @param max_search_frames Maximum frames to search in each direction
 * @param prefer_direction Preferred search direction: -1=backward, 0=either, 1=forward
 * @return Frame position of nearest zero crossing, or target_frame if none found
 */
static sf_count_t find_nearest_zero_crossing(const AudioData *audio, sf_count_t target_frame,
                                            sf_count_t max_search_frames, int prefer_direction)
{
    if (!audio || !audio->samples || target_frame >= audio->frame_count) {
        return target_frame;
    }
    
    sf_count_t search_start = (target_frame > max_search_frames) ? 
                              target_frame - max_search_frames : 0;
    sf_count_t search_end = (target_frame + max_search_frames < audio->frame_count) ? 
                            target_frame + max_search_frames : audio->frame_count - 1;
    
    sf_count_t best_crossing = target_frame;
    sf_count_t best_distance = max_search_frames + 1;
    
    /* Search for zero crossings */
    for (sf_count_t i = search_start; i < search_end; i++) {
        float current_sample = 0.0f;
        float next_sample = 0.0f;
        
        /* Get current and next sample (mono conversion if needed) */
        if (audio->channels == 2) {
            float left = audio->samples[i * 2];
            float right = audio->samples[i * 2 + 1];
            current_sample = (left + right) / 2.0f;
            
            if (i + 1 < audio->frame_count) {
                left = audio->samples[(i + 1) * 2];
                right = audio->samples[(i + 1) * 2 + 1];
                next_sample = (left + right) / 2.0f;
            }
        } else {
            current_sample = audio->samples[i];
            if (i + 1 < audio->frame_count) {
                next_sample = audio->samples[i + 1];
            }
        }
        
        /* Check for zero crossing between current and next sample */
        if (fabsf(current_sample) < MIN_ZERO_CROSSING_AMPLITUDE ||
            (current_sample * next_sample <= 0.0f && 
             fabsf(current_sample) + fabsf(next_sample) > MIN_ZERO_CROSSING_AMPLITUDE)) {
            
            sf_count_t distance = (i > target_frame) ? i - target_frame : target_frame - i;
            
            /* Apply direction preference */
            int matches_preference = 1;
            if (prefer_direction < 0 && i > target_frame) {
                matches_preference = 0;  /* Prefer backward but this is forward */
            } else if (prefer_direction > 0 && i < target_frame) {
                matches_preference = 0;  /* Prefer forward but this is backward */
            }
            
            /* Update best crossing if this is better */
            if (matches_preference && distance < best_distance) {
                best_crossing = i;
                best_distance = distance;
            } else if (!matches_preference && distance < best_distance && best_distance > max_search_frames) {
                /* Only use non-preferred direction if no preferred direction found */
                best_crossing = i;
                best_distance = distance;
            }
        }
    }
    
    return best_crossing;
}

/**
 * Find the attack point (rapid energy increase) near a target frame
 *
 * @param audio Audio data to search in
 * @param target_frame Target frame position (beat detection point)
 * @param search_window_frames Window size to search for attack
 * @return Frame position of attack onset, or target_frame if none found
 */
static sf_count_t find_attack_point(const AudioData *audio, sf_count_t target_frame, 
                                   sf_count_t search_window_frames)
{
    if (!audio || !audio->samples || target_frame >= audio->frame_count) {
        return target_frame;
    }
    
    /* Define energy calculation window (small, ~2ms) */
    sf_count_t energy_window = (sf_count_t)(0.002 * audio->sample_rate); /* 2ms window */
    if (energy_window < 32) energy_window = 32; /* Minimum 32 samples */
    
    sf_count_t search_start = (target_frame > search_window_frames) ? 
                              target_frame - search_window_frames : 0;
    sf_count_t search_end = (target_frame + search_window_frames < audio->frame_count) ? 
                            target_frame + search_window_frames : audio->frame_count - 1;
    
    sf_count_t best_attack_frame = target_frame;
    double max_energy_ratio = 1.0;
    
    /* Scan through the search window */
    for (sf_count_t i = search_start; i < search_end - energy_window; i++) {
        /* Calculate energy before and after this point */
        double energy_before = 0.0;
        double energy_after = 0.0;
        
        /* Energy before (preceding window) */
        for (sf_count_t j = 0; j < energy_window && i - j > 0; j++) {
            sf_count_t frame_idx = i - j;
            float sample = 0.0f;
            
            if (audio->channels == 2) {
                float left = audio->samples[frame_idx * 2];
                float right = audio->samples[frame_idx * 2 + 1];
                sample = (left + right) / 2.0f;
            } else {
                sample = audio->samples[frame_idx];
            }
            
            energy_before += sample * sample;
        }
        
        /* Energy after (following window) */
        for (sf_count_t j = 0; j < energy_window && i + j < audio->frame_count; j++) {
            sf_count_t frame_idx = i + j;
            float sample = 0.0f;
            
            if (audio->channels == 2) {
                float left = audio->samples[frame_idx * 2];
                float right = audio->samples[frame_idx * 2 + 1];
                sample = (left + right) / 2.0f;
            } else {
                sample = audio->samples[frame_idx];
            }
            
            energy_after += sample * sample;
        }
        
        /* Calculate energy ratio */
        double energy_ratio = (energy_before > 0.0001) ? energy_after / energy_before : 1.0;
        
        /* Check if this represents a significant attack */
        if (energy_ratio > ATTACK_ENERGY_THRESHOLD && energy_ratio > max_energy_ratio) {
            max_energy_ratio = energy_ratio;
            best_attack_frame = i;
        }
    }
    
    return best_attack_frame;
}

/**
 * Find optimal slice boundary combining attack detection and zero-crossing alignment
 *
 * @param audio Audio data to search in
 * @param beat_frame Beat detection frame position
 * @param is_start_boundary True if this is a slice start, false if slice end
 * @return Optimal frame position for slice boundary
 */
static sf_count_t find_optimal_slice_boundary(const AudioData *audio, sf_count_t beat_frame, 
                                             int is_start_boundary)
{
    if (!audio || !audio->samples) {
        return beat_frame;
    }
    
    /* Calculate search windows */
    sf_count_t attack_window = (sf_count_t)(ATTACK_DETECTION_WINDOW_MS * audio->sample_rate / 1000.0);
    sf_count_t zero_search_window = (sf_count_t)(MAX_ZERO_CROSSING_SEARCH_MS * audio->sample_rate / 1000.0);
    sf_count_t max_pre_attack = (sf_count_t)(MAX_PRE_ATTACK_MS * audio->sample_rate / 1000.0);
    
    if (is_start_boundary) {
        /* For slice starts: Find attack point, then nearby zero crossing */
        
        /* Step 1: Find the attack point */
        sf_count_t attack_frame = find_attack_point(audio, beat_frame, attack_window);
        
        /* Step 2: Find zero crossing near attack point, but constrained by pre-attack limit */
        sf_count_t min_allowed_frame = (attack_frame > max_pre_attack) ? 
                                      attack_frame - max_pre_attack : 0;
        
        sf_count_t zero_frame = find_nearest_zero_crossing(audio, attack_frame, 
                                                          zero_search_window, -1); /* Prefer backward */
        
        /* Step 3: Ensure we don't go too far before the attack */
        if (zero_frame < min_allowed_frame) {
            zero_frame = min_allowed_frame;
        }
        
        /* Step 4: Verify we're not cutting off the attack */
        if (zero_frame > attack_frame) {
            zero_frame = attack_frame; /* Don't cut into the attack */
        }
        
        return zero_frame;
        
    } else {
        /* For slice ends: Use standard zero crossing with forward preference */
        return find_nearest_zero_crossing(audio, beat_frame, zero_search_window, 1);
    }
}

int load_audio_file(const char *filename, AudioData *audio_data) 
{
    SNDFILE *sf_file;
    SF_INFO sf_info;
    
    if (!filename || !audio_data) {
        LOG_ERROR("Invalid parameters to load_audio_file");
        return BEATSLICE_ERROR_INVALID_PARAMETER;
    }
    
    /* Initialize SF_INFO structure */
    memset(&sf_info, 0, sizeof(SF_INFO));
    
    /* Open the audio file */
    sf_file = sf_open(filename, SFM_READ, &sf_info);
    if (!sf_file) {
        LOG_ERROR("Failed to open audio file '%s': %s", filename, sf_strerror(NULL));
        return BEATSLICE_ERROR_FILE_NOT_FOUND;
    }
    
    /* Validate audio format */
    if (sf_info.frames <= 0 || sf_info.samplerate <= 0 || sf_info.channels <= 0) {
        LOG_ERROR("Invalid audio format in file '%s'", filename);
        sf_close(sf_file);
        return BEATSLICE_ERROR_UNSUPPORTED_FORMAT;
    }
    
    if (sf_info.channels > 2) {
        LOG_WARNING("Multi-channel audio detected (%d channels), will be mixed to mono/stereo", 
                   sf_info.channels);
    }
    
    /* Allocate memory for audio samples */
    audio_data->frame_count = sf_info.frames;
    audio_data->sample_rate = sf_info.samplerate;
    audio_data->channels = sf_info.channels;
    audio_data->duration = (double)sf_info.frames / sf_info.samplerate;
    audio_data->sf_info = sf_info;
    
    size_t total_samples = sf_info.frames * sf_info.channels;
    audio_data->samples = malloc(sizeof(float) * total_samples);
    
    if (!audio_data->samples) {
        LOG_ERROR("Failed to allocate memory for audio samples (%zu bytes)", 
                 sizeof(float) * total_samples);
        sf_close(sf_file);
        return BEATSLICE_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Read audio samples */
    sf_count_t frames_read = sf_readf_float(sf_file, audio_data->samples, sf_info.frames);
    
    if (frames_read != sf_info.frames) {
        LOG_WARNING("Expected %ld frames, read %ld frames", 
                   (long)sf_info.frames, (long)frames_read);
        audio_data->frame_count = frames_read;
        audio_data->duration = (double)frames_read / sf_info.samplerate;
    }
    
    sf_close(sf_file);
    
    LOG_INFO("Loaded audio: %.2fs, %d Hz, %d channel(s), %ld frames",
             audio_data->duration, audio_data->sample_rate, 
             audio_data->channels, (long)audio_data->frame_count);
    
    return BEATSLICE_SUCCESS;
}

int export_beat_slices(const AudioData *audio, const AnalysisResults *analysis,
                      const BeatsliceConfig *config) 
{
    if (!audio || !analysis || !config) {
        LOG_ERROR("Invalid parameters to export_beat_slices");
        return BEATSLICE_ERROR_INVALID_PARAMETER;
    }
    
    if (analysis->beat_count <= 0) {
        LOG_WARNING("No beats detected, creating single slice of entire audio");
    }
    
    int sf_format, sf_subtype;
    if (get_output_format(config->output_format, &sf_format, &sf_subtype) != 0) {
        LOG_ERROR("Unsupported output format: %s", config->output_format);
        return BEATSLICE_ERROR_INVALID_PARAMETER;
    }
    
    int slices_exported = 0;
    
    /* Calculate attack detection parameters */
    sf_count_t attack_window_frames = (sf_count_t)(ATTACK_DETECTION_WINDOW_MS * audio->sample_rate / 1000.0);
    
    LOG_DEBUG("Attack-focused alignment: attack_window=%.1fms, max_pre_attack=%.1fms", 
              ATTACK_DETECTION_WINDOW_MS, MAX_PRE_ATTACK_MS);
    
    /* Calculate number of slices: beats + 1 (or 1 if no beats) */
    int slice_count = (analysis->beat_count > 0) ? analysis->beat_count + 1 : 1;
    
    /* Create slice information array */
    SliceInfo *slices = malloc(sizeof(SliceInfo) * slice_count);
    if (!slices) {
        LOG_ERROR("Failed to allocate slice information array");
        return BEATSLICE_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Calculate slice boundaries with attack-focused alignment */
    for (int i = 0; i < slice_count; i++) {
        SliceInfo *slice = &slices[i];
        double original_start_time, original_end_time;
        sf_count_t original_start_frame, original_end_frame;
        sf_count_t aligned_start_frame, aligned_end_frame;
        
        if (analysis->beat_count == 0) {
            /* No beats detected - single slice of entire audio */
            original_start_time = 0.0;
            original_end_time = audio->duration;
            slice->confidence = 0.0;
        } else if (i == 0) {
            /* First slice: from start to first beat */
            original_start_time = 0.0;
            original_end_time = analysis->beat_times[0];
            slice->confidence = analysis->confidence[0];
        } else if (i == analysis->beat_count) {
            /* Last slice: from last beat to end */
            original_start_time = analysis->beat_times[analysis->beat_count - 1];
            original_end_time = audio->duration;
            slice->confidence = analysis->confidence[analysis->beat_count - 1];
        } else {
            /* Middle slices: between consecutive beats */
            original_start_time = analysis->beat_times[i - 1];
            original_end_time = analysis->beat_times[i];
            slice->confidence = analysis->confidence[i - 1];
        }
        
        /* Convert original times to frame indices */
        original_start_frame = (sf_count_t)(original_start_time * audio->sample_rate);
        original_end_frame = (sf_count_t)(original_end_time * audio->sample_rate);
        
        /* Find attack-optimized boundaries */
        aligned_start_frame = find_optimal_slice_boundary(audio, original_start_frame, 1); /* Start boundary */
        aligned_end_frame = find_optimal_slice_boundary(audio, original_end_frame, 0);     /* End boundary */
        
        /* Set final slice boundaries */
        slice->start_frame = aligned_start_frame;
        slice->frame_count = aligned_end_frame - aligned_start_frame;
        slice->start_time = (double)aligned_start_frame / audio->sample_rate;
        slice->end_time = (double)aligned_end_frame / audio->sample_rate;
        
        /* Calculate boundary adjustments for debug output */
        double start_adjustment = slice->start_time - original_start_time;
        double end_adjustment = slice->end_time - original_end_time;
        
        /* Calculate pre-attack duration for start boundaries */
        double pre_attack_ms = 0.0;
        if (i > 0) { /* Not the first slice */
            sf_count_t attack_frame = find_attack_point(audio, original_start_frame, attack_window_frames);
            double attack_time = (double)attack_frame / audio->sample_rate;
            pre_attack_ms = (attack_time - slice->start_time) * 1000.0;
        }
        
        /* Validate slice length - preserve boundary slices */
        double slice_duration = slice->end_time - slice->start_time;
        int is_boundary_slice = (i == 0 || i == slice_count - 1);
        
        if (!is_boundary_slice && slice_duration < config->min_slice_length) {
            LOG_DEBUG("Skipping middle slice %d (%.3fs < %.3fs minimum)", 
                     i, slice_duration, config->min_slice_length);
            slice->frame_count = 0;  /* Mark as invalid */
            continue;
        }
        
        /* Always preserve boundary slices, even if they're very short */
        if (is_boundary_slice && slice_duration < config->min_slice_length) {
            LOG_DEBUG("Preserving boundary slice %d despite short duration (%.3fs)", 
                     i, slice_duration);
        }
        
        /* Ensure minimum frame count for valid slices */
        if (slice->frame_count <= 0) {
            LOG_WARNING("Slice %d has no frames after attack optimization, skipping", i);
            continue;
        }
        
        /* Generate output filename */
        snprintf(slice->filename, sizeof(slice->filename), "%s_%03d.%s",
                config->output_prefix, i, config->output_format);
        
        if (config->verbose) {
            LOG_DEBUG("Slice %d: %.3fs - %.3fs (%.3fs, %ld frames, conf: %.2f)", 
                     i, slice->start_time, slice->end_time, slice_duration,
                     (long)slice->frame_count, slice->confidence);
            
            /* Show attack-focused adjustments */
            if (fabs(start_adjustment) > 0.001 || fabs(end_adjustment) > 0.001) {
                if (i > 0 && pre_attack_ms >= 0.0) {
                    LOG_DEBUG("  Attack-focused alignment: start %+.1fms, end %+.1fms, pre-attack %.1fms", 
                             start_adjustment * 1000.0, end_adjustment * 1000.0, pre_attack_ms);
                } else {
                    LOG_DEBUG("  Boundary alignment: start %+.1fms, end %+.1fms", 
                             start_adjustment * 1000.0, end_adjustment * 1000.0);
                }
            }
            
            /* Warn if pre-attack time is excessive */
            if (i > 0 && pre_attack_ms > MAX_PRE_ATTACK_MS) {
                LOG_DEBUG("  WARNING: Pre-attack time (%.1fms) exceeds maximum (%.1fms) - may have silence before hit", 
                         pre_attack_ms, MAX_PRE_ATTACK_MS);
            }
        }
    }
    
    /* Export each valid slice */
    for (int i = 0; i < slice_count; i++) {
        SliceInfo *slice = &slices[i];
        
        if (slice->frame_count <= 0) {
            continue;  /* Skip invalid slices */
        }
        
        /* Ensure output directory exists */
        if (ensure_output_directory(slice->filename) != 0) {
            LOG_ERROR("Failed to create output directory for %s", slice->filename);
            continue;
        }
        
        /* Set up output file info */
        SF_INFO out_info = audio->sf_info;
        out_info.format = sf_format | sf_subtype;
        out_info.frames = slice->frame_count;
        
        /* Open output file */
        SNDFILE *out_file = sf_open(slice->filename, SFM_WRITE, &out_info);
        if (!out_file) {
            LOG_ERROR("Failed to create output file '%s': %s", 
                     slice->filename, sf_strerror(NULL));
            continue;
        }
        
        /* Calculate source data pointer and validate bounds */
        sf_count_t max_start_frame = audio->frame_count - 1;
        if (slice->start_frame > max_start_frame) {
            LOG_ERROR("Invalid start frame %ld for slice %d (max: %ld)", 
                     (long)slice->start_frame, i, (long)max_start_frame);
            sf_close(out_file);
            continue;
        }
        
        /* Adjust frame count if it would exceed audio bounds */
        sf_count_t available_frames = audio->frame_count - slice->start_frame;
        if (slice->frame_count > available_frames) {
            LOG_WARNING("Adjusting slice %d frame count from %ld to %ld (audio boundary)", 
                       i, (long)slice->frame_count, (long)available_frames);
            slice->frame_count = available_frames;
            out_info.frames = slice->frame_count;
        }
        
        float *source_data = audio->samples + (slice->start_frame * audio->channels);
        
        /* Copy and optionally normalize slice data */
        if (config->normalize_output) {
            /* Allocate temporary buffer for normalization */
            size_t slice_samples = slice->frame_count * audio->channels;
            float *slice_buffer = malloc(sizeof(float) * slice_samples);
            
            if (slice_buffer) {
                memcpy(slice_buffer, source_data, sizeof(float) * slice_samples);
                normalize_audio_slice(slice_buffer, slice->frame_count, audio->channels);
                
                sf_count_t frames_written = sf_writef_float(out_file, slice_buffer, slice->frame_count);
                free(slice_buffer);
                
                if (frames_written != slice->frame_count) {
                    LOG_WARNING("Expected to write %ld frames, wrote %ld frames to %s",
                               (long)slice->frame_count, (long)frames_written, slice->filename);
                }
            } else {
                LOG_WARNING("Failed to allocate normalization buffer, writing unnormalized");
                sf_writef_float(out_file, source_data, slice->frame_count);
            }
        } else {
            /* Write slice data directly */
            sf_count_t frames_written = sf_writef_float(out_file, source_data, slice->frame_count);
            
            if (frames_written != slice->frame_count) {
                LOG_WARNING("Expected to write %ld frames, wrote %ld frames to %s",
                           (long)slice->frame_count, (long)frames_written, slice->filename);
            }
        }
        
        sf_close(out_file);
        slices_exported++;
        
        if (config->verbose) {
            LOG_INFO("Exported attack-optimized slice %d: %s (%.3fs, confidence: %.2f)",
                    i, slice->filename, slice->end_time - slice->start_time, slice->confidence);
        }
    }
    
    free(slices);
    
    LOG_INFO("Exported %d attack-optimized audio slices", slices_exported);
    return slices_exported;
}

int export_analysis_data(const AnalysisResults *analysis, const BeatsliceConfig *config) 
{
    if (!analysis || !config) {
        return BEATSLICE_ERROR_INVALID_PARAMETER;
    }
    
    char json_filename[512];
    snprintf(json_filename, sizeof(json_filename), "%s_analysis.json", config->output_prefix);
    
    FILE *json_file = fopen(json_filename, "w");
    if (!json_file) {
        LOG_ERROR("Failed to create analysis export file '%s': %s", 
                 json_filename, strerror(errno));
        return BEATSLICE_ERROR_EXPORT_FAILED;
    }
    
    /* Write JSON analysis data */
    fprintf(json_file, "{\n");
    fprintf(json_file, "  \"version\": \"%d.%d.%d\",\n", 
            BEATSLICE_VERSION_MAJOR, BEATSLICE_VERSION_MINOR, BEATSLICE_VERSION_PATCH);
    fprintf(json_file, "  \"analysis\": {\n");
    fprintf(json_file, "    \"beat_count\": %d,\n", analysis->beat_count);
    fprintf(json_file, "    \"estimated_bpm\": %.2f,\n", analysis->estimated_bpm);
    fprintf(json_file, "    \"tempo_confidence\": %.3f,\n", analysis->tempo_confidence);
    fprintf(json_file, "    \"frame_count\": %d,\n", analysis->frame_count);
    fprintf(json_file, "    \"hop_time\": %.6f,\n", analysis->analysis_hop_time);
    fprintf(json_file, "    \"attack_optimized_slicing\": true\n");
    fprintf(json_file, "  },\n");
    
    /* Export beat times */
    fprintf(json_file, "  \"beats\": [\n");
    for (int i = 0; i < analysis->beat_count; i++) {
        fprintf(json_file, "    {\"time\": %.6f, \"confidence\": %.3f}%s\n",
                analysis->beat_times[i], analysis->confidence[i],
                (i < analysis->beat_count - 1) ? "," : "");
    }
    fprintf(json_file, "  ],\n");
    
    /* Export energy envelope (sampled for size) */
    fprintf(json_file, "  \"energy_envelope\": [");
    int envelope_step = fmax(1, analysis->frame_count / 1000);  /* Max 1000 points */
    for (int i = 0; i < analysis->frame_count; i += envelope_step) {
        if (i > 0) fprintf(json_file, ", ");
        fprintf(json_file, "%.4f", analysis->energy_envelope ? analysis->energy_envelope[i] : 0.0);
    }
    fprintf(json_file, "],\n");
    
    /* Export onset function (sampled) */
    fprintf(json_file, "  \"onset_function\": [");
    for (int i = 0; i < analysis->frame_count; i += envelope_step) {
        if (i > 0) fprintf(json_file, ", ");
        fprintf(json_file, "%.4f", analysis->onset_function ? analysis->onset_function[i] : 0.0);
    }
    fprintf(json_file, "]\n");
    
    fprintf(json_file, "}\n");
    
    fclose(json_file);
    
    LOG_INFO("Exported analysis data to %s", json_filename);
    return BEATSLICE_SUCCESS;
}

const char* get_output_filename(const char *prefix, int slice_index, const char *format) 
{
    static char filename[512];
    snprintf(filename, sizeof(filename), "%s_%03d.%s", prefix, slice_index, format);
    return filename;
}

/*----- Private Functions -----*/

static int get_output_format(const char *format_name, int *sf_format, int *sf_subtype) 
{
    for (int i = 0; output_formats[i].extension != NULL; i++) {
        if (strcmp(format_name, output_formats[i].extension) == 0) {
            *sf_format = output_formats[i].sf_format;
            *sf_subtype = output_formats[i].sf_subtype;
            return 0;
        }
    }
    return -1;
}

static int ensure_output_directory(const char *filepath) 
{
    char *filepath_copy = strdup(filepath);
    char *dir_path = dirname(filepath_copy);
    
    struct stat st;
    int result = 0;
    
    /* Check if directory exists */
    if (stat(dir_path, &st) == 0) {
        if (!S_ISDIR(st.st_mode)) {
            LOG_ERROR("Output path exists but is not a directory: %s", dir_path);
            result = -1;
        }
    } else {
        /* Create directory */
        if (mkdir(dir_path, 0755) != 0 && errno != EEXIST) {
            LOG_ERROR("Failed to create output directory '%s': %s", 
                     dir_path, strerror(errno));
            result = -1;
        }
    }
    
    free(filepath_copy);
    return result;
}

static void normalize_audio_slice(float *samples, sf_count_t frame_count, int channels) 
{
    if (!samples || frame_count <= 0) {
        return;
    }
    
    /* Find peak amplitude */
    float peak = 0.0f;
    sf_count_t total_samples = frame_count * channels;
    
    for (sf_count_t i = 0; i < total_samples; i++) {
        float abs_sample = fabsf(samples[i]);
        if (abs_sample > peak) {
            peak = abs_sample;
        }
    }
    
    /* Normalize to 0.95 to avoid clipping */
    if (peak > 0.0f) {
        float scale = 0.95f / peak;
        for (sf_count_t i = 0; i < total_samples; i++) {
            samples[i] *= scale;
        }
    }
}

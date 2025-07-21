/*
 * fileio.h - Audio file I/O operations header
 * 
 * Function declarations and constants for audio file loading and exporting
 * using libsndfile with support for multiple formats.
 */

#ifndef FILEIO_H
#define FILEIO_H

#include <sndfile.h>
#include "beatslice.h"

/* File format constants */
#define MAX_FILENAME_LENGTH     512
#define MAX_FORMAT_NAME_LENGTH  16
#define MAX_PATH_LENGTH         1024

/* Supported output format information */
typedef struct {
    const char *name;           /* Format name (e.g., "wav", "aiff") */
    const char *description;    /* Human-readable description */
    int sf_format;             /* libsndfile format constant */
    int sf_subtype;            /* libsndfile subtype constant */
    const char *extension;     /* File extension */
    int max_channels;          /* Maximum supported channels */
    int supports_float;        /* Whether format supports floating point */
} AudioFormatInfo;

/* File validation results */
typedef struct {
    int is_valid;              /* Whether file is valid audio */
    int is_supported;          /* Whether format is supported */
    double duration_seconds;   /* File duration in seconds */
    int64_t file_size_bytes;   /* File size in bytes */
    char format_name[64];      /* Detected format name */
    char error_message[256];   /* Error description if invalid */
} FileValidationResult;

/* Export statistics */
typedef struct {
    int slices_attempted;      /* Number of slices attempted to export */
    int slices_successful;     /* Number of slices successfully exported */
    int slices_skipped;        /* Number of slices skipped (too short, etc.) */
    double total_duration;     /* Total duration of exported audio */
    int64_t total_bytes;       /* Total bytes written */
    double processing_time;    /* Time taken for export operation */
} ExportStatistics;

/*----- Core Audio I/O Functions -----*/

/**
 * Load an audio file into memory
 *
 * Loads an audio file using libsndfile, supporting multiple formats.
 * The audio data is converted to floating-point and stored in the
 * provided AudioData structure.
 *
 * @param filename Path to the audio file to load
 * @param audio_data Pointer to AudioData structure to populate
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int load_audio_file(const char *filename, AudioData *audio_data);

/**
 * Export beat slices as separate audio files
 *
 * Takes the original audio data and analysis results to create
 * individual audio files for each detected beat segment.
 *
 * @param audio Pointer to original audio data
 * @param analysis Pointer to analysis results with beat positions
 * @param config Pointer to configuration (output format, etc.)
 * @return Number of slices exported, or negative error code
 */
int export_beat_slices(const AudioData *audio, const AnalysisResults *analysis,
                      const BeatsliceConfig *config);

/**
 * Export analysis data to JSON file
 *
 * Exports beat detection analysis data to a JSON file for visualization
 * or further processing.
 *
 * @param analysis Pointer to analysis results
 * @param config Pointer to configuration (output prefix, etc.)
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int export_analysis_data(const AnalysisResults *analysis, const BeatsliceConfig *config);

/**
 * Get formatted output filename for a slice
 *
 * Generates a properly formatted filename for an audio slice with
 * the specified prefix, index, and format.
 *
 * @param prefix Output filename prefix
 * @param slice_index Index of the slice (zero-based)
 * @param format Output format extension
 * @return Pointer to static filename string
 */
const char* get_output_filename(const char *prefix, int slice_index, 
                               const char *format);

/*----- Advanced File Operations -----*/

/**
 * Validate audio file before processing
 *
 * Checks if a file exists, is readable, and contains valid audio data
 * that can be processed by the beat detection system.
 *
 * @param filename Path to the audio file
 * @param result Pointer to FileValidationResult structure
 * @return BEATSLICE_SUCCESS if validation completed (check result.is_valid)
 */
int validate_audio_file(const char *filename, FileValidationResult *result);

/**
 * Get supported audio format information
 *
 * Returns information about supported audio formats for input and output.
 *
 * @param format_name Format name to query (e.g., "wav", "flac")
 * @return Pointer to AudioFormatInfo structure, or NULL if not found
 */
const AudioFormatInfo* get_format_info(const char *format_name);

/**
 * List all supported audio formats
 *
 * Returns an array of all supported audio format information.
 *
 * @param count Pointer to store the number of formats
 * @return Array of AudioFormatInfo structures
 */
const AudioFormatInfo* get_supported_formats(int *count);

/**
 * Export single audio slice
 *
 * Exports a single audio slice with specified parameters. Used internally
 * by export_beat_slices but can be called directly for custom slicing.
 *
 * @param audio Source audio data
 * @param start_time Start time in seconds
 * @param end_time End time in seconds
 * @param filename Output filename
 * @param format Output format info
 * @param normalize Whether to normalize the output
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int export_single_slice(const AudioData *audio, double start_time, double end_time,
                       const char *filename, const AudioFormatInfo *format,
                       int normalize);

/**
 * Create audio preview/overview file
 *
 * Creates a compressed audio file with beat markers for preview purposes.
 * Useful for visualization and verification of beat detection results.
 *
 * @param audio Original audio data
 * @param analysis Analysis results with beat positions
 * @param filename Output filename for preview
 * @param duration Maximum duration for preview (0 = full length)
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int create_audio_preview(const AudioData *audio, const AnalysisResults *analysis,
                        const char *filename, double duration);

/*----- File Utility Functions -----*/

/**
 * Check if output directory exists and is writable
 *
 * Verifies that the directory for output files exists and can be written to.
 * Creates the directory if it doesn't exist.
 *
 * @param filepath Full path to a file (directory will be extracted)
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int ensure_output_directory_exists(const char *filepath);

/**
 * Get audio file metadata
 *
 * Extracts metadata from an audio file without loading the entire file.
 *
 * @param filename Path to audio file
 * @param sample_rate Pointer to store sample rate
 * @param channels Pointer to store channel count
 * @param frames Pointer to store frame count
 * @param duration Pointer to store duration in seconds
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int get_audio_metadata(const char *filename, int *sample_rate, int *channels,
                      sf_count_t *frames, double *duration);

/**
 * Check available disk space for output
 *
 * Estimates required disk space for export operation and checks if
 * sufficient space is available.
 *
 * @param output_dir Output directory path
 * @param estimated_size Estimated space needed in bytes
 * @return 1 if sufficient space available, 0 otherwise
 */
int check_available_disk_space(const char *output_dir, int64_t estimated_size);

/**
 * Generate unique filename if file exists
 *
 * If the specified filename already exists, generates a unique variant
 * by appending a number (e.g., file_001.wav, file_002.wav).
 *
 * @param base_filename Base filename to use
 * @param unique_filename Buffer to store unique filename
 * @param buffer_size Size of unique_filename buffer
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int generate_unique_filename(const char *base_filename, char *unique_filename,
                            size_t buffer_size);

/*----- Audio Format Conversion Functions -----*/

/**
 * Convert audio sample format
 *
 * Converts audio samples between different bit depths and formats.
 *
 * @param input_samples Source audio samples
 * @param output_samples Destination buffer
 * @param sample_count Number of samples to convert
 * @param input_format Source format (SF_FORMAT_*)
 * @param output_format Destination format (SF_FORMAT_*)
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int convert_sample_format(const void *input_samples, void *output_samples,
                         sf_count_t sample_count, int input_format, int output_format);

/**
 * Resample audio data
 *
 * Changes the sample rate of audio data using high-quality resampling.
 * Note: This function requires additional resampling library.
 *
 * @param input_audio Source audio data
 * @param output_audio Destination audio data structure
 * @param target_sample_rate Desired output sample rate
 * @return BEATSLICE_SUCCESS on success, error code on failure
 */
int resample_audio(const AudioData *input_audio, AudioData *output_audio,
                  int target_sample_rate);

/*----- Error Handling and Diagnostics -----*/

/**
 * Get detailed libsndfile error information
 *
 * Returns a detailed error description from libsndfile for the last
 * operation that failed.
 *
 * @param sf_file SNDFILE handle (can be NULL for global errors)
 * @return String describing the last error
 */
const char* get_sndfile_error_string(SNDFILE *sf_file);

/**
 * Print supported formats information
 *
 * Prints a formatted list of all supported audio formats to stdout.
 * Useful for help messages and debugging.
 */
void print_supported_formats(void);

/**
 * Get file size in bytes
 *
 * Returns the size of a file in bytes, or -1 if the file doesn't exist
 * or cannot be accessed.
 *
 * @param filename Path to the file
 * @return File size in bytes, or -1 on error
 */
int64_t get_file_size(const char *filename);

/**
 * Check if filename has valid audio extension
 *
 * Checks if a filename has an extension that corresponds to a supported
 * audio format.
 *
 * @param filename Filename to check
 * @return 1 if valid audio extension, 0 otherwise
 */
int has_valid_audio_extension(const char *filename);

/*----- Constants for Format Support -----*/

/* Maximum supported values */
#define MAX_SUPPORTED_SAMPLE_RATE  192000
#define MIN_SUPPORTED_SAMPLE_RATE  8000
#define MAX_SUPPORTED_CHANNELS     8
#define MAX_SUPPORTED_BIT_DEPTH    32

/* File size limits (in bytes) */
#define MAX_INPUT_FILE_SIZE        (2LL * 1024 * 1024 * 1024)  /* 2GB */
#define MIN_INPUT_FILE_SIZE        1024                        /* 1KB */

/* Export quality settings */
#define DEFAULT_EXPORT_QUALITY     0.8    /* Compression quality 0.0-1.0 */
#define DEFAULT_NORMALIZATION_PEAK 0.95   /* Peak level for normalization */

#endif /* FILEIO_H */

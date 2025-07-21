# Beatslice Algorithm: Detailed Technical Explanation

## Overview

Beatslice uses a **spectral energy-based onset detection** algorithm that analyzes frequency domain characteristics to identify transient events (beats) in audio. The core principle is that beat onsets typically cause sudden increases in spectral energy, particularly in the low-frequency range where drum kicks and bass elements reside.

## Phase 1: Initialization and Configuration

### Default Parameters
```c
#define DEFAULT_FFT_SIZE       2048    // 2048 samples = ~46ms at 44.1kHz
#define DEFAULT_HOP_SIZE       512     // 512 samples = ~11.6ms at 44.1kHz
#define DEFAULT_SENSITIVITY    0.5     // Mid-range sensitivity
#define DEFAULT_MIN_BPM        60.0    // 60 BPM = 1 beat per second
#define DEFAULT_MAX_BPM        200.0   // 200 BPM = 3.33 beats per second
```

### Key Configuration Elements

- **FFT Size (2048)**: Determines frequency resolution and time window
  - Frequency resolution = sample_rate / fft_size = 44100/2048 = ~21.5 Hz per bin
  - Time window = fft_size / sample_rate = 2048/44100 = ~46ms

- **Hop Size (512)**: Controls overlap and temporal resolution
  - 75% overlap (512/2048 = 0.25, so 75% overlap)
  - Analysis frame rate = sample_rate / hop_size = 44100/512 = ~86 frames/second
  - Time between analyses = ~11.6ms

## Phase 2: Audio Loading and Preprocessing

### Sample Format Handling
```c
// CRITICAL: Maintains libsndfile's normalized float range (-1.0 to +1.0)
// No scaling to 16-bit integer range - keeps modern float precision
```

### Stereo to Mono Conversion
```c
if (audio->channels == 2) {
    // Average left and right channels
    double left = audio->samples[(offset + i) * 2];
    double right = audio->samples[(offset + i) * 2 + 1];
    output[i] = (left + right) / 2.0;
}
```

## Phase 3: Spectral Analysis Pipeline

### 3.1 Windowing Process

For each analysis frame:

1. **Extract Audio Window**: Extract `fft_size` samples starting at `frame_index * hop_size`

2. **Apply Hann Window**: Reduces spectral leakage
```c
// Hann window: w[n] = 0.5 * (1 - cos(2π*n/(N-1)))
multiplier = 0.5 * (1.0 - cos(2.0 * M_PI * i / (size - 1)));
buffer[i] *= multiplier;
```

3. **FFT Computation**: Uses FFTW3 for efficient real-to-complex FFT
```c
fftw_execute(g_analysis_state.fft_plan);
// Output: fft_output[0] to fft_output[fft_size/2] (1025 bins for 2048 FFT)
```

### 3.2 Frequency Band Analysis

The algorithm divides the spectrum into three frequency bands:

#### Low Frequency Band (Bass/Kick Detection)
- **Range**: 40 Hz - 250 Hz
- **Purpose**: Captures kick drums, bass drops, low-frequency transients
- **Bin Range**: bins 2-12 (at 44.1kHz, 2048 FFT)
- **Primary detection mechanism**

#### Mid Frequency Band
- **Range**: 250 Hz - 2000 Hz  
- **Purpose**: Snares, vocals, mid-range instruments
- **Bin Range**: bins 12-93

#### High Frequency Band
- **Range**: 2000 Hz - Nyquist (22.05 kHz)
- **Purpose**: Hi-hats, cymbals, high-frequency content
- **Bin Range**: bins 93-1024

### 3.3 Energy Calculation

For each frequency band:
```c
// Calculate magnitude from complex FFT output
real = fft_out[i][0];  // Real component
imag = fft_out[i][1];  // Imaginary component
magnitude = sqrt(real*real + imag*imag);
energy += magnitude;

// Normalize by number of bins
energy /= (high_bin - low_bin + 1);
```

### 3.4 Additional Features (Optional)

#### Zero Crossing Rate (if enabled)
```c
// Counts sign changes in the time domain
if (prev_sign != 0 && curr_sign != 0 && prev_sign != curr_sign) {
    crossings++;
}
zero_crossing_rate = crossings / window_size;
```

#### RMS Energy
```c
// Root Mean Square energy calculation
sum += sample * sample;
rms_energy = sqrt(sum / window_size);
```

## Phase 4: Beat Detection Algorithm

### 4.1 Energy Envelope Creation

- Uses **low-frequency energy** as the primary onset detection function
- Normalizes envelope to 0-1 range for consistent thresholding
- Stores energy value for each analysis frame (~86 frames/second)

### 4.2 Adaptive Thresholding (if enabled)

Creates a dynamic threshold that adapts to local signal characteristics:

```c
// Sliding window average with scaling factor
window_size = min_peak_distance * 2;  // Typically ~25 frames
threshold_factor = 0.3 + (sensitivity * 0.5);  // Range: 0.3 to 0.8

for each frame:
    local_average = average_energy_in_window(frame, window_size);
    adaptive_threshold[frame] = local_average * threshold_factor;
```

### 4.3 Sensitivity Mapping

User sensitivity (0.0-1.0) maps to algorithm parameters:

#### Fixed Threshold Mode
```c
threshold = 0.15 + (sensitivity * 0.25);
// sensitivity 0.0 → threshold 0.15 (less sensitive)
// sensitivity 1.0 → threshold 0.40 (more sensitive)
```

#### Adaptive Threshold Mode
```c
peak_threshold_factor = 0.3 + (sensitivity * 0.5);
// sensitivity 0.0 → factor 0.30
// sensitivity 1.0 → factor 0.80
```

### 4.4 Peak Detection Algorithm

#### Minimum Peak Distance
```c
min_peak_distance_sec = 0.3;  // 300ms minimum between beats
min_peak_distance_frames = min_peak_distance_sec * sample_rate / hop_size;
// At 44.1kHz: 0.3 * 44100 / 512 ≈ 26 frames
```

#### Peak Finding Logic
```c
for each frame i in energy_envelope:
    if (envelope[i] > threshold &&           // Exceeds threshold
        envelope[i] > envelope[i-1] &&       // Rising edge
        envelope[i] >= envelope[i+1]) {      // Peak or plateau
        
        // Check minimum distance constraint
        valid = true;
        for each existing_peak:
            if (abs(i - existing_peak) < min_peak_distance) {
                if (envelope[i] > envelope[existing_peak]) {
                    replace existing_peak with i;  // Keep higher peak
                } else {
                    valid = false;  // Skip this peak
                }
            }
        
        if (valid) {
            add i to peak_list;
        }
    }
```

## Phase 5: Beat Time Conversion

### Frame Index to Time Conversion
```c
beat_times[i] = peak_frame_index * hop_size / sample_rate;
// Example: frame 86 at 44.1kHz, hop 512
// time = 86 * 512 / 44100 = 1.0 seconds
```

### BPM Estimation
```c
if (beat_count > 1) {
    total_time = last_beat_time - first_beat_time;
    estimated_bpm = (beat_count - 1) * 60.0 / total_time;
}
```

## Phase 6: Slice Boundary Creation

### Slice Definition
- **First slice**: 0.0 → first_beat_time
- **Middle slices**: beat[i-1] → beat[i]  
- **Last slice**: last_beat_time → end_of_audio

### Validation
- Minimum slice length filter (default 100ms)
- Boundary slices (first/last) preserved regardless of length
- Frame boundary validation

## Algorithm Characteristics

### Temporal Resolution
- **Analysis rate**: ~86 Hz (every 11.6ms)
- **Minimum beat spacing**: 300ms (200 BPM max)
- **Typical accuracy**: ±11.6ms (one hop size)

### Frequency Resolution  
- **Bin width**: ~21.5 Hz
- **Low-freq focus**: 2-12 bins (40-250 Hz)
- **Total spectrum**: 1025 bins (0-22.05 kHz)

### Latency Characteristics
- **Look-ahead**: 46ms (FFT window size)
- **Processing delay**: Minimal (real-time capable)
- **Total system latency**: ~57ms (window + hop)

### Strengths
- Excellent for percussive, electronic music
- Robust low-frequency onset detection
- Adaptive to varying signal levels
- Computationally efficient

### Limitations
- May miss soft onsets without strong spectral change
- Low-frequency focus can miss high-frequency-only onsets
- Fixed frequency band boundaries
- Requires sufficient SNR in low frequencies

## Tuning Parameters

### For Different Music Styles

**Electronic/Dance Music** (default settings work well):
- FFT: 2048, Hop: 512, Sensitivity: 0.5-0.7

**Acoustic/Jazz** (softer onsets):
- FFT: 1024, Hop: 256, Sensitivity: 0.3-0.5
- Enable zero crossings, use adaptive threshold

**Hip-Hop/Drums** (strong transients):
- FFT: 2048, Hop: 512, Sensitivity: 0.6-0.8
- Higher threshold values

**Classical** (complex harmonic content):
- FFT: 4096, Hop: 1024, Sensitivity: 0.2-0.4
- Longer analysis windows for stability

This algorithm represents a well-balanced approach to onset detection, combining proven DSP techniques with practical real-world performance considerations.
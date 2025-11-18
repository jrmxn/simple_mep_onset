# MEP Onset Detection (MATLAB)

This repository provides a compact MATLAB implementation for detecting motor evoked potential (MEP) onset latency from EMG-like waveforms. The main entry point is:

```matlab
t_onset_ms = get_valid_onset( ...
    y_mep, t, latency_est_auc, t_onset_bounds, fs, ...
    th_onset_sd, wiggle_percentage_cutoff, th_onset_uv);
```
It returns the latency in milliseconds or `NaN` if a valid onset is not found.

---

## Quick start

1. Place `get_valid_onset.m` in your MATLAB path.
2. Call it with a vector signal `y_mep` and a matching time vector `t` (in seconds).

Minimal example with synthetic data:
```matlab
fs = 5000;                 % Hz
dt = 1/fs;
t  = (-0.05:dt:0.20)';     % seconds
rng(1);
y  = 5e-6*randn(size(t));  % ~5 uV noise

onset_true = 0.028;        % seconds
y = y + (t >= onset_true) .* (50e-6 * exp(-((t - onset_true)/0.010)));

latency_est_auc = NaN;     % kept for API compatibility, not used
t_onset_bounds  = [0, 0.060];    % seconds; samples before t_onset_bounds(1) are blanked
th_onset_sd     = 4;       % multiplier on baseline SD
wiggle_percentage_cutoff = 2;  % percent of global |y|
th_onset_uv     = 20e-6;   % 20 uV gate for the first valid sample

t_onset_ms = get_valid_onset( ...
    y, t, latency_est_auc, t_onset_bounds, fs, ...
    th_onset_sd, wiggle_percentage_cutoff, th_onset_uv);

fprintf('Detected onset: %.2f ms\n', t_onset_ms);
```

---

## Function reference

### `get_valid_onset(y_mep, t, latency_est_auc, t_onset_bounds, fs, th_onset_sd, wiggle_percentage_cutoff, th_onset_uv)`

Detects the earliest sustained excursion above a baseline-derived threshold and returns its onset time in **milliseconds**. A set of validity checks must pass or the function returns `NaN`.

| Argument | Type | Units | Description |
|---|---|---|---|
| `y_mep` | vector | signal units (e.g., uV) | EMG or MEP waveform. |
| `t` | vector | seconds | Time vector, same length as `y_mep`. |
| `latency_est_auc` | scalar | n/a | **Not used** in this implementation. Kept for API compatibility. |
| `t_onset_bounds` | 1x2 vector | seconds | `[t_min, t_max]`. All samples with `t < t_min` are hard-blanked to 0 before detection. |
| `fs` | scalar | Hz | Sampling rate. |
| `th_onset_sd` | scalar | multiplier | Threshold is `th_onset_sd * sd_baseline`, where `sd_baseline` is computed from `|y_mep|` for `t < -1 ms`. |
| `wiggle_percentage_cutoff` | scalar | percent | If the initial deflection after threshold crossing is smaller than this percent of the global `|y|`, the detector skips forward to the next deflection to avoid tiny wiggles. Use `NaN` or `Inf` to disable. |
| `th_onset_uv` | scalar | signal units | Absolute amplitude gate that the first nonzero post-blanking sample must be below. Guards against large artifacts. |

**Return**: `t_onset_ms` (scalar). Onset latency in milliseconds, or `NaN` if no valid onset is found.

**Conventions and guards**
- Baseline standard deviation is estimated from `t < -1e-3` seconds on `abs(y_mep)` with `omitnan`.
- Excursion must remain above threshold for about **2 ms** (controlled internally).
- Onset must occur earlier than **50 ms**. Later detections are treated as likely artifacts.
- If indices needed for validation are missing, the function returns `NaN`.

---

## How it works

1. **Blanking**: All samples where `t < t_onset_bounds(1)` are set to zero. This suppresses the stimulation artifact window.
2. **Baseline SD**: Compute `sd_baseline = std(abs(y_mep(t < -1 ms)), 'omitnan')`.
3. **Stable excursion**: Find the earliest time where `|y|` exceeds `th_onset_sd * sd_baseline` and remains above this threshold for ~2 ms.
4. **Traceback**: Smooth the envelope and trace back to the start of the deflection to refine the onset index.
5. **Wiggle skip**: Optionally skip tiny pre-onset wiggles based on a percentage of the global maximum absolute amplitude.
6. **Validity checks**: Enforce early-onset gate (< 50 ms), small first nonzero sample, and other simple guards. Otherwise return `NaN`.

The helper routines are defined in the same `.m` file and are considered implementation details:
- `get_onset_from_stable_excursion`
- `skip_wiggle`
- `onset_traceback`
- `filtfilt_smooth`

---

## Parameter tips

- `t_onset_bounds`: start with `[0, 0.06]` seconds for TMS MEPs. Increase the upper bound if you expect very late responses.
- `th_onset_sd`: 3 to 6 is a reasonable range depending on SNR.
- `wiggle_percentage_cutoff`: 1 to 5 typically suppresses tiny pre-onset bumps without hiding real onsets.
- `th_onset_uv`: pick a value near the noise floor of your preprocessing chain, for example 10 to 30 uV if your signal is in microvolts.

---

## Assumptions and limitations

- `t` is aligned so that `t = 0` corresponds to the stimulus time. If not, align before calling.
- The detector expects a single dominant deflection that rises above the threshold and remains there briefly.
- Baseline noise should be stationary in the pre-stimulus window used for SD estimation.
- The code uses simple moving average smoothing for the traceback step and does not remove NaNs. Mirror the preprocessing you use in analysis for best results.

---

## Troubleshooting

- **Always `NaN`**: Check that `t` and `y_mep` have the same length, `fs` is correct, and `t_onset_bounds(1)` is not negative relative to your time base. Verify that you have nonzero data after blanking.
- **Late artifacts flagged as onset**: Reduce `max_latency_ms` inside the file or raise `th_onset_sd` and `th_onset_uv`.
- **Missed weak onsets**: Lower `th_onset_sd` and consider reducing `wiggle_percentage_cutoff`.
- **Multiple small bumps**: Increase `wiggle_percentage_cutoff` or the internal stable excursion window (currently ~2 ms).

---

## Citation

Was used in:  https://doi.org/10.1113/JP286183

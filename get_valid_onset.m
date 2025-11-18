function t_onset_th = get_valid_onset( ...
    y_mep, t, latency_est_auc, t_onset_bounds, fs, ...
    th_onset_sd, wiggle_percentage_cutoff, th_onset_uv)
%GET_VALID_ONSET Detect MEP onset (ms) with simple validity checks.
%
% Parameters
%   y_mep : vector, EMG/MEP signal (arbitrary units, e.g., uV)
%   t     : vector, time in seconds (same length as y_mep)
%   latency_est_auc : unused (kept for API compatibility)
%   t_onset_bounds  : [t_min, t_max] seconds; samples with t < t_min are blanked
%   fs              : sampling rate (Hz)
%   th_onset_sd     : threshold multiplier on baseline SD (abs signal)
%   wiggle_percentage_cutoff : percent of global |y| to ignore tiny pre-onset wiggles
%   th_onset_uv     : absolute amplitude threshold at first valid sample (units of y_mep)
%
% Returns
%   t_onset_th : onset latency in ms (NaN if invalid/undetected)
%
% Notes
%   - Baseline SD is computed over t < -1 ms (abs signal).
%   - Onset must occur < 50 ms to be considered valid.
%   - Excursion must remain over threshold for ~2 ms.


% ---- constants
t_onset_th = nan;
t_stable_excursion = 2e-3;         % seconds of stable excursion
baseline_limit_s   = -1e-3;        % seconds (baseline window upper bound)
max_latency_ms     = 50;           % ms (late artifact cutoff)
first_nonzero_uv   = 2.0;          % absolute amplitude gate for first nonzero

% ---- basic blanking before earliest allowed onset
pre_stim = t < t_onset_bounds(1);
y_mep_blanked         = y_mep;
y_mep_blanked(pre_stim) = 0;
y_abs_blanked         = abs(y_mep_blanked);

% First valid (non-zero) sample index after blanking
ix_first_valid = find(y_abs_blanked, 1, 'first');

% ---- baseline SD (use abs signal, pre -1 ms, omit NaNs)
sd_baseline = std(abs(y_mep(t < baseline_limit_s)), 'omitnan');

% ---- thresholded onset with stable excursion
ix_onset_th = get_onset_from_stable_excursion( ...
    y_mep_blanked, sd_baseline * th_onset_sd, t_stable_excursion * fs, ...
    fs, wiggle_percentage_cutoff);

% ---- validity checks (structured to avoid empty-index errors)
is_valid_estimate = true;

% 1) when data first appears after stim, it should be below an ABS threshold
if ~isempty(ix_first_valid)
    is_valid_estimate = is_valid_estimate & (y_abs_blanked(ix_first_valid) < th_onset_uv);
else
    is_valid_estimate = false;
end

% 2) onset must be early enough (late artifacts rejected)
if ~isempty(ix_onset_th)
    is_valid_estimate = is_valid_estimate & (t(ix_onset_th) * 1e3 < max_latency_ms);
else
    is_valid_estimate = false;
end

% 3) the very first nonzero sample in the (blanked) record should be small
ix_first_nonzero = find(y_mep_blanked, 1, 'first');
if ~isempty(ix_first_nonzero)
    is_valid_estimate = is_valid_estimate & (abs(y_mep_blanked(ix_first_nonzero)) < first_nonzero_uv);
else
    is_valid_estimate = false;
end

% ---- output (ms)
if ~isempty(ix_onset_th) && is_valid_estimate
    t_onset_th = t(ix_onset_th) * 1e3;
end
end


% ======================================================================
% Helpers
% ======================================================================
function ix_onset = get_onset_from_stable_excursion(y, threshold, n_stable_excursion, fs, wiggle_percentage_cutoff)
%GET_ONSET_FROM_STABLE_EXCURSION
%   Returns the index where |y| first exceeds 'threshold' and remains
%   above it for ~n_stable_excursion samples, with traceback and wiggle-skip.

y_abs = abs(y);
case_excursion = (y_abs > threshold);

% Candidate points where excursion is true
vec_onset = find(case_excursion);

% Require excursion to last for n_stable_excursion samples
% (moving average of the boolean mask near 1 means all were true)
vec_excursion_lasting = filter(ones(1, n_stable_excursion) / n_stable_excursion, 1, case_excursion);
vec_onset_lasting    = find(vec_excursion_lasting > 1 - eps);

if isempty(vec_onset_lasting)
    ix_onset = [];
    return;
end

% Find the first onset that leads into a lasting region
ix_onset = [];
for k = 1:numel(vec_onset)
    cand = vec_onset(k);
    ix_lasting = vec_onset_lasting(vec_onset_lasting > cand);
    if isempty(ix_lasting)
        continue; % no lasting region after this candidate
    end
    ix_lasting = ix_lasting(1);
    if all(vec_excursion_lasting(cand:ix_lasting) > 0)
        ix_onset = cand;
        break;
    end
end

if isempty(ix_onset)
    return;
end

% Trace back to the start of this deflection on a smoothed envelope
ix_onset = onset_traceback(y_abs, ix_onset, fs);

% Optionally skip tiny initial wiggles near the threshold
if ~isempty(ix_onset)
    max_wiggle_skip = 5;
    n_skips = 0;
    while true
        ix_new = skip_wiggle(y, ix_onset, wiggle_percentage_cutoff);
        if isempty(ix_new) || ix_new == ix_onset || n_skips > max_wiggle_skip
            break;
        end
        ix_onset = ix_new;
        n_skips = n_skips + 1;
    end
end

% Guard: onset must not be too close to the end of the record
if ~isempty(ix_onset) && ix_onset >= (length(y_abs) - n_stable_excursion)
    ix_onset = [];
end
end


function ix_onset = skip_wiggle(y, ix_onset, wiggle_percentage_cutoff)
%SKIP_WIGGLE Skip small pre-onset wiggles based on sign changes and size.

if ~isfinite(wiggle_percentage_cutoff) || isempty(ix_onset)
    return; % no change
end

dy_sign = diff(sign(y));
% mask everything up to (and including) the current onset
clip_to = min(numel(dy_sign), ix_onset + 1);
dy_sign(1:clip_to) = NaN;

ix_second_onset_th = find(abs(dy_sign) > 0, 1, 'first');
if isempty(ix_second_onset_th)
    return; % nothing to skip
end

% size of first deflection vs global max
max_initial_deflect = max(abs(y(ix_onset:ix_second_onset_th)));
global_max          = max(abs(y));
if global_max == 0
    return;
end
proportion_of_max = max_initial_deflect / global_max;

if proportion_of_max < (wiggle_percentage_cutoff * 1e-2)
    ix_onset = ix_second_onset_th;
end
end


function ix_onset_early = onset_traceback(y, ix_onset, fs)
%ONSET_TRACEBACK Move onset backward to the start of the deflection.
n = max(1, round(fs * 0.5e-3));           % ~0.5 ms smoothing
yf = filtfilt_smooth(y, n);               % zero-phase moving average

d = sign([0, diff(yf)]);                  % slope sign
d(ix_onset:end) = 0;

% detect last transition into a rising segment before ix_onset
m  = fliplr(filter([+1, -1], 1, fliplr(d)));
mk = (m == -2);
ix_onset_early = find(mk, 1, 'last');
end


function yf = filtfilt_smooth(y, n)
%FILTFILT_SMOOTH Zero-phase moving-average smoothing (handles NaNs like original).
b  = ones(n, 1) / n;
yf = fliplr(filter(b, 1, fliplr(filter(b, 1, y))));
end

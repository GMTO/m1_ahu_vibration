%
%
%

clearvars
% Save PSD into a MATLAB data file
save_psd = false;
% Remove frequency points outside the definied range
rm_freq_pts = true;


%% Load data
%%
input_dt_fname = 'fan_wrench_reconstructed_psd_timeseries_60sec.mat';
load(input_dt_fname,'timeseriesAHUfm', 'time', 'psdAHUf', 'freq')


%% Calculate the PSD
%%
% Sampling parameters
dt = time(2) - time(1);
fs = 1 / dt;
N = length(time);

% Target a fine frequency resolution (~0.1 Hz) to isolate motor harmonics
target_df = 0.1;
nfft = 2^nextpow2(fs / target_df);
nfft = min(nfft, 2^floor(log2(N))); % Bound by total signal length

% Use Hann window for correct energy/variance integration
window = hann(nfft);
noverlap = floor(nfft * 0.5); % 50% overlap

% Pre-allocate 3D matrix
S_uu = zeros(6, 6, floor(nfft/2)+1);

% Populating ALL 6x6 cross and auto-spectral terms
for i = 1:6
    for j = 1:6
        [S_uu(i,j,:), freq_s_uu] = cpsd(timeseriesAHUfm(:,i), timeseriesAHUfm(:,j), ...
            window, noverlap, nfft, fs);
    end
end

%% PSD Comparison (main diagonal terms)
%%
ch_names = {'F_x', 'F_y', 'F_z', 'M_x', 'M_y', 'M_z'};

figure(101);
set(gcf,"Position",[831   667   919   590], 'Name', 'PSD Comparison')
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for ch = 1:6
    nexttile
    loglog(freq, psdAHUf(:, ch), 'b', 'LineWidth', 1.4); hold on;
    loglog(freq_s_uu, squeeze(S_uu(ch, ch, :)), '--', 'LineWidth', 1.5);
    xlim([min(freq), max(freq)]); grid on; hold off;

    title(sprintf('%s', ch_names{ch}));

    % Adjust units based on force vs. moment
    if ch <= 3, ylabel('PSD (N^2/Hz)');
    else, ylabel('PSD ((N\cdotm)^2/Hz)'); xlabel('Frequency (Hz)');
    end
    % Legend
    if ch == 1, legend('Original psdAHUf', 'Recomputed (Welch)'); end
end

% Constrain frequency range
if(rm_freq_pts)
    freq_mask = (freq_s_uu >= 10) & (freq_s_uu <= 200);
    freq_s_uu = freq_s_uu(freq_mask);
    S_uu = S_uu(:, :, freq_mask);
end


%% Check if the PSD matrix is Hermitian positive semi-definite
%%
% Audit S_uu for non-Hermitian or negative eigenvalue issues
neg_eig_count = 0;
for f = 1:length(freq_s_uu)
    % Force symmetry
    S_f = (S_uu(:,:,f) + S_uu(:,:,f)') / 2;
    e = eig(S_f);
    if any(e < -1e-12), neg_eig_count = neg_eig_count + 1; end
end
fprintf('\n-> Frequencies with negative eigenvalues in S_uu: %d / %d.\n',...
    neg_eig_count, length(freq_s_uu));

if(neg_eig_count > 0)
    warning('Regularizing input load PSD');
    % Clean and regularize S_uu matrix across all frequency points
    for f = 1:length(freq_s_uu)
        % Force exact Hermitian symmetry
        S_f = (S_uu(:,:,f) + S_uu(:,:,f)') / 2;
        [V, D] = eig(S_f);
        d = diag(D);

        % Only rebuild if negative numerical noise exists
        if any(d < 0)
            d(d < 0) = 0;
            S_f = V * diag(d) * V';
        end

        S_uu(:,:,f) = S_f;
    end
end


%% Save full vibration input load PSD
%%
if(save_psd)
    psd_fname = 'fan_wrench_rec_psd_20250826_rro';
    save(psd_fname,'input_dt_fname','freq_s_uu','S_uu');
    fprintf("\n-> Input load PSD data saved in \n%s.mat\n",...
        fullfile(pwd,psd_fname));
end


%% Time- & Frequency-domain standard deviation comparison.
% Direct Time-Domain Standard Deviation
ch = 1;
std_time = std(timeseriesAHUfm(:, ch));

% Frequency-Domain Integrated Standard Deviation
% Extract auto-PSD for Channel ch from S_uu
S_11 = squeeze(S_uu(ch, ch, :));
df = freq_s_uu(2) - freq_s_uu(1);
std_freq = sqrt(trapz(freq_s_uu, S_11));

% Calculate Input Load Ratio
ratio_input = std_time / std_freq;

fprintf('\n --- Time- & Frequency-domain standard deviation comparison ---\n')
fprintf('Time-domain Force Std:       %.4f\n', std_time);
fprintf('Freq-domain Integrated Std:  %.4f\n', std_freq);
fprintf('Input Load Mismatch Ratio:   %.4f\n', ratio_input);


%% Coherence Analysis (Zoomed 10–200 Hz)
%%
try close(102); catch, end
figure(102);
set(gcf,"Position",[431   267   919   590],...
    'Name', 'Motor Vibration Load Cross-Coupling (10-200 Hz)')
coherence_matrix = zeros(floor(nfft/2)+1, 6, 6);
for i = 1:6
    for j = i+1:6
        [Cij, f_coh] = mscohere(timeseriesAHUfm(:, i), timeseriesAHUfm(:, j), window, noverlap, nfft, fs);
        coherence_matrix(:, i, j) = Cij;
        coherence_matrix(:, j, i) = Cij;
        
        semilogx(f_coh, Cij, 'DisplayName', sprintf('%s \\leftrightarrow %s', ch_names{i}, ch_names{j}));
        if(i==1 && j==2), hold on; end
    end
end


grid on; hold off;
xlim([10, 200]); %ylim([0, 1]);
xlabel('Frequency (Hz)');
ylabel('Coherence \gamma^2(f)');
title('Cross-Coupling Coherence Spectrum (10–200 Hz)');
legend('Location', 'northeastoutside');

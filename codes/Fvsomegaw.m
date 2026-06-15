%{
Fvsomegaw.m

This code creates a plot showing how the size of the faithful set decreases
with increasing omegaw.

This code uses the Communication Toolbox.

Last modified by David A. Edwards on 6/14/26.

%}

% Define parameters.

fmax = 512;                    % desired analog bandwidth (Hz)
% i: index of rows
n = 512; % Size of Hadamard matrix
% omegaw: Word transmission rate.
omegawmax = 64; % Maximum word transmission rate.

% Define vectors to be used.

bestvec = zeros(n,1); % Vector of row of best matches.
faithful = zeros(n,1); % Vector of faithful set.
plotset = zeros(omegawmax,2); % Plotting vector

% Digital -> Analog -> Digital demo
% rng(1);
% n = 512;
% x = 2*(rand(1,n) > 0.5) - 1;   % bipolar ±1 digital symbols

% omegaw = 1;

for ffactor = 1:2

    fmax = 512*ffactor

for omegaw = 1:omegawmax

    for i = 1:n
        % See if this can be vectorized!
        % Compute the best match for row i.
        bestvec(i) = bestmatch(i,n,omegaw,fmax);
    end

    % Echo value
    if mod(omegaw,8)==0 
        fprintf('omegaw = %d\n', omegaw);
    end

    % Compute the faithful set.
    faithful = (bestvec == (1:length(bestvec))');
    % Compute the size of the faithful set, and store for later plotting:
    plotset(omegaw,ffactor)=sum(faithful);

end

end

% Then plot the result:
figure;
hold on;

% Column 1: black circles connected
h1 = plot(plotset(:,1), '-o', 'Color', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

% Column 2: red triangles connected
h2 = plot(plotset(:,2), '-^', 'Color', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

title('|F| vs. word rate');
xlabel('Word rate');
ylabel('|F|');
xlim([1 omegawmax]);
ylim([0 n]); 


% Create legend labels "fmax = 512*column number"
labels = arrayfun(@(c) sprintf('Analog cutoff = %d', 512*c), 1:2, 'UniformOutput', false);
legend([h1 h2], labels, 'Location', 'best');
grid on;

hold off;



% omegawmax



function f = bestmatch(i,n,omegaw,fmax)
% This function computes the index of the closest transmitted encoding to the received
% encoding i.

% Called by: main
% Calls: none.

% Input parameters:
% fmax: desired analog bandwidth (Hz)
% i: input row number
% n: size of Hadamard matrix
% omegaw: word transmission rate.

% Output variables:
% bestmatch: index of encoding which best matches the received state

% Internal variables:
fs_sym = n*omegaw;                 % symbol/sample rate (Hz) - must be > 2*fmax_symbol? here choose >= 2*fmax/n/A
fs_analog = fs_sym;              % analog sampling frequency (Hz), must satisfy fs_analog > 2*fmax
h = hadamard(n); % Hadamard matrix
x = h(i,:); % encoding to be transmitted.

% Upsample ratio p/q for resample: convert from fs_sym to fs_analog
% p = fs_analog;
% q = fs_sym;

% Create analog (bandlimited) waveform by resampling (interpolation with antialias filter)
% resample returns a sequence sampled at fs_analog
analog = resample(x,fs_analog,fs_sym);    % anti-aliasing/interpolation filter applied

% Time vectors
t_sym = (0:n-1)/fs_sym;
t_analog = (0:length(analog)-1)/fs_analog;

% Verify analog bandwidth (optional): design lowpass to enforce fmax if needed
% Here resample's built-in filter already limits to nyquist of symbol rate, but to ensure fmax use filtfilt:
Wn = fmax/(fs_analog/2);       % normalized cutoff for analog sampling
if Wn < 1
    [b,a] = butter(6, Wn);     % 6th-order Butterworth lowpass
    analog = filtfilt(b,a,double(analog));
end

% Recover digital by sampling analog at symbol instants (nearest indices)
L = fs_analog / fs_sym;        % integer upsample factor (should be integer)
if abs(L - round(L)) > 1e-10
    error('fs_analog must be an integer multiple of fs_sym for simple downsampling. Use resample for arbitrary ratios.');
end
L = round(L);
recovered_samples = analog(1:L:end);   % pick samples corresponding to symbol instants

% Decision device (hard decision to ±1)
x_rec = sign(recovered_samples);
x_rec(x_rec==0) = 1;           % tie-break if exact zero

% Now that we have computed x_rec, we find the encoding with minimum
% distance from it.
d = sum(h ~= x_rec, 2);    % m-by-1 vector of Hamming distances from each row
% Here the second returned output is the index of the minimum, which is
% what we need.
[min_value, f] = min(d);

% Metrics
% mse_analog_vs_sym = mean((interp1(t_sym, x, t_analog, 'nearest') - analog).^2); % analog vs nearest-step
% mse_rec = mean((x - x_rec).^2);
% bit_errors = sum(x ~= x_rec);
% ber = bit_errors / n;

% fprintf('MSE (after analog interp vs nearest step) = %.6f\n', mse_analog_vs_sym);
% fprintf('MSE (recovered symbols) = %.6f\n', mse_rec);
%  fprintf('Bit errors = %d, BER = %.6f\n', bit_errors, ber);

end

% % Plots: show a short segment for clarity
% Lseg = min(100, n);
% idx_sym = 1:Lseg;
% idx_analog = 1:round(L*Lseg);
% 
% figure;
% subplot(3,1,1);
% stem(t_sym(idx_sym), x(idx_sym), 'b', 'filled');
% title('Original Digital Symbols (±1)');
% xlabel('Time (s)'); ylabel('Amplitude');
% xlim([t_sym(1) t_sym(Lseg)]); ylim([-1.5 1.5]); grid on;
% 
% subplot(3,1,2);
% plot(t_analog(1:idx_analog(end)), analog(1:idx_analog(end)), 'k-');
% hold on;
% stairs(t_sym(idx_sym), x(idx_sym), 'b--','LineWidth',1);
% title('Analog Waveform (bandlimited interpolation) and Original Symbols');
% xlabel('Time (s)'); ylabel('Amplitude');
% xlim([t_sym(1) t_sym(Lseg)]); ylim([-1.5 1.5]); legend('Analog','Symbols'); grid on;
% 
% subplot(3,1,3);
% stem(t_sym(idx_sym), x_rec(idx_sym), 'r', 'filled');
% hold on;
% stem(t_sym(idx_sym), x(idx_sym), 'b','filled');
% title('Recovered Digital Symbols vs Original (short segment)');
% xlabel('Time (s)'); ylabel('Amplitude');
% xlim([t_sym(1) t_sym(Lseg)]); ylim([-1.5 1.5]); legend('Recovered','Original'); grid on;

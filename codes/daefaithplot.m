%{
daefaithplot.m

This code performs several faithful-encoding calculations to provide the
related figures in the report.

This code uses the Signal Processing Toolbox, 
Communication Toolbox and the Curve Fitting Toolbox.

Last modified by David A. Edwards on 6/16/26.

%}


% Define parameters.

fmax = 512; % desired analog bandwidth (Hz): also called omegaa.
% i: index of rows
n = 512; % Size of Hadamard matrix
% omegaw: Word transmission rate.
omegawmax = 64; % Maximum word transmission rate.

% Variables
bestvec = zeros(n,1); % Vector of row of best matches.
% faithful: Vector of faithful set.
plotset = zeros(omegawmax,2); % Plotting vector

% Number of runs to use with noise.
maxrun = 50;
maxrun = 10;

% Define vectors to be used.

% omegaw = 1;

% % Set LaTeX interpreters and default font size for all plots.
set(groot, 'defaultTextInterpreter', 'latex', ...
           'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex', ...
           'defaultAxesFontSize', 14, ...
           'defaultTextFontSize', 14); 

%

% First plot: F vs omegaw, no noise.

% snr = Inf corresponds to no noise.
snr = Inf;

for ffactor = 1:2

    fmax = n*ffactor

    for omegaw = 1:omegawmax

        faithful = getfaithful(n,omegaw,fmax,snr);

        % For a particular value of omegaw, keep faithful for quoting
        if omegaw==30
            f30 = faithful;
        end
        % Compute the size of the faithful set, and store for later plotting.
        % We want to do an exponential fit, so we take the log.
        plotset(omegaw,ffactor)=log(sum(faithful));

    end

end

% % Then plot the result:

omegaplot(plotset,snr,n);

%

% Second plot: add noise, but just for n=512.

% Define plot arrays.
% faithmat contains the results from all the runs.  We use a cell array so
% that we can store the arrays inside without worrying about three indices:
faithmat = cell(maxrun,1);
% noisyplot contains the columns to be plotted.  We only use the fourth
% column in the other plot.
noisyplot = zeros(omegawmax,4);

% Loop over snr values.

for snr = -3:-6:-9

    for omegaw = 1:omegawmax

        for run = 1:maxrun

            faithful = getfaithful(n,omegaw,fmax,snr);
            faithmat{run,1} = faithful;

        end

        % Compute the size of the faithful set, and store for later plotting.
        % We want to do an exponential fit, so we take the log.
        % The first entry is one run picked at random.
        randrun = randi(maxrun);
        noisyplot(omegaw,1)=log(sum(faithmat{randrun,1}));
        % The second entry is the average.  First we sum each of the runs to
        % find F:
        test = cellfun(@sum,faithmat);
        % Then take the mean of all those sums:
        noisyplot(omegaw,2)=log(mean(test));
        % The last entry is the intersection.  So we count the number of rows
        % where EVERY entry is 1:

        % faithmat is maxrun-by-1 cell, each cell is a 1-by-N or N-by-1 vector of 0/1
        M = vertcat(faithmat{:});        % maxrun-by-N numeric matrix

        colsAllOnes = all(M == 1, 1);    % 1-by-N logical: true where every run has a 1
       % scalar: number of indices i where all cells have 1
        noisyplot(omegaw,3)=log(sum(colsAllOnes));

    end

    % % Then plot the result:

    nomegaplot(noisyplot,snr,n,maxrun);

end

%%
% Third plot: F vs bandwidth, no noise.

omegaw = 10;
% This counts the row number because fmax jumps.
count = 1;
snr = Inf;
for fmax = 512:8:1024 % Go from 512 to 1024 by 8s

    faithful = getfaithful(n,omegaw,fmax,snr);

    % Echo value
    if mod(fmax,64)==0
        fprintf('fmax = %d\n', fmax);
    end

    % Compute the size of the faithful set, and store for later plotting:
    plotset(count,:)=[fmax,sum(faithful)];
    count = count+1;

end

% % Then plot the result:

bandplot(plotset,snr,n,omegaw);

% Fourth plot: F vs bandwidth, with noise.
for snr = -3:-6:-9
% for snr = -3:-3
    % Reset the row number.
    count = 1;
    for fmax = 512:8:1024 % Go from 512 to 1024 by 8s

        for run = 1:maxrun

            faithful = getfaithful(n,omegaw,fmax,snr);
            faithmat{run,1} = faithful;

        end

        % Echo value
        if mod(fmax,64)==0
            fprintf('fmax = %d\n', fmax);
        end

        % Compute the size of the faithful set, and store for later plotting.
        % We want to do an exponential fit, so we take the log.
        % The first entry is one run picked at random.
        randrun = randi(maxrun);
        noisyplot(count,1)=sum(faithmat{randrun,1});
        % The second entry is the average.  First we sum each of the runs to
        % find F:
        test = cellfun(@sum,faithmat);
        % Then take the mean of all those sums:
        noisyplot(count,2)=mean(test);
        % The last entry is the intersection.  So we count the number of rows
        % where EVERY entry is 1:

        % faithmat is maxrun-by-1 cell, each cell is a 1-by-N or N-by-1 vector of 0/1
        M = vertcat(faithmat{:});        % maxrun-by-N numeric matrix

        colsAllOnes = all(M == 1, 1);    % 1-by-N logical: true where every run has a 1
        % scalar: number of indices i where all cells have 1

        noisyplot(count,3)=sum(colsAllOnes);

        % Compute the size of the faithful set, and store for later plotting:
        noisyplot(count,4)=fmax;
        count = count+1;

    end

    % % Then plot the result:

    nbandplot(noisyplot,snr,n,maxrun,omegaw);

end


%%
function f = bestmatch(i,h,n,omegaw,fmax,snr)
% This function computes the index of the closest transmitted encoding to the
% received
% encoding i.

% Called by: faithful
% Calls: adam_DAD, dae_DAD.

% Input parameters:
% fmax: desired analog bandwidth (Hz)
% i: input row number
% h: Walsh matrix
% n: size of Walsh matrix
% omegaw: word transmission rate.
% snr: signal-to-noise ratio (dB)

% Output variables:
% bestmatch: index of encoding which best matches the received state

% Internal variables:
% DADflag: 1 if using Matlab's converter; 2 if using Adam's

DADflag = 1;

% Internal variables:
fs_sym = n*omegaw;  % symbol/sample rate (Hz) - must be > 2*fmax_symbol? here choose >= 2*fmax/n/A
fs_analog = fs_sym; % analog sampling frequency (Hz), must satisfy fs_analog > 2*fmax
% Create the nxn Hadamard matrix in Walsh order


% Upsample ratio p/q for resample: convert from fs_sym to fs_analog
% p = fs_analog;
% q = fs_sym;

% The input string is the ith row of the Walsh matrix:
x = h(i,:);

% Now compute the received signal using one of the two possible methods:
if DADflag==1
    % Copilot method
    x_rec = dae_DAD(x,fs_analog,fs_sym,n,fmax,snr);
else
    % Adam's method.  IMPORTANT: The second argument isn't right.
    x_rec = adam_DAD(x',fs_analog/omegaw^2);
end

% x_rec;

% Now that we have computed x_rec, we find the encoding with minimum
% distance from it.
d = sum(h ~= x_rec, 2);    % m-by-1 vector of Hamming distances from each row
% Next we find the indices of all entries that match the minimum:
f = find(d == min(d));
% For the encoding to be faithful, it must be the ONLY minimum, so we check
% how many indices are listed.  So if we have more than one minimum, we
% just set it to zero so it doesn't match:
if length(f)~=1
    f = 0;
end

% Metrics
% mse_analog_vs_sym = mean((interp1(t_sym, x, t_analog, 'nearest') - analog).^2); % analog vs nearest-step
% mse_rec = mean((x - x_rec).^2);
% bit_errors = sum(x ~= x_rec);
% ber = bit_errors / n;

% fprintf('MSE (after analog interp vs nearest step) = %.6f\n', mse_analog_vs_sym);
% fprintf('MSE (recovered symbols) = %.6f\n', mse_rec);
%  fprintf('Bit errors = %d, BER = %.6f\n', bit_errors, ber);

end

function return_data = adam_DAD(inp_seq,cutoff)

% This function takes an input sequence and given a bandwidth cutoff,
% returns an output sequence.  This version was written by Adam Petrucci
% using only the FFT.  This version can't handle noise.

% Called by: bestmatch
% Calls: none

% Input parameters:
% cutoff: bandwidth cutoff (omega_a)
% inp_seq: input sequence

% Output parameters:
% return_data: sequence fed through DAD process

% choose code length
code_length = length(inp_seq);

% rudimentary discrete ADPCM (only at intervals)
inp_con = [0;cumsum(inp_seq)];

% fill in intervals
m = 64;
inp_con_fill = repelem(inp_con,m);
M = length(inp_con_fill);

% visualize cont input
%inp_anal_fig = figure(1);
%clf(inp_anal_fig);
%inp_anal_ax = axes(inp_anal_fig);
%plot(inp_anal_ax,1:M,inp_con_fill);
%title(inp_anal_ax,'Input Analog');

% Define frequencies and cutoff
k = [0:M/2-1 -M/2:-1]';
xi = 2*pi*k/code_length;

% fft -> filter -> ifft
Fhat = fft(inp_con_fill);
mask = abs(xi) <= cutoff;
Fhat_filt = Fhat .* mask;
out_con_fill = real(ifft(Fhat_filt));

% visualize cont output
%out_anal_fig = figure(2);
%clf(out_anal_fig);
%out_anal_ax = axes(out_anal_fig);
%plot(out_anal_ax,1:M,out_con_fill);
%title(out_anal_ax,'Output Analog');

% extract values at intervals
out_cont_fill = out_con_fill(1:m:end);
out_seq = 2*(diff(out_cont_fill) >= 0) - 1;

return_data = out_seq;

%inp_seq.'
%out_seq.'

%num_diff = sum(inp_seq ~= out_seq);

%return_data = num_diff;

end

function x_rec = dae_DAD(x,fs_analog,fs_sym,n,fmax,snr)

% This function takes and given parameters about the low-pass cutoff,
% returns an output sequence.  This version was written by David A. Edwards
% and Copilot and uses several "black-box" Maple functions.

% Called by: bestmatch
% Calls: none

% Input parameters:
% fmax: desired analog bandwidth (Hz)
% fs_sym = n*omegaw: symbol/sample rate (Hz) - must be > 2*fmax_symbol? here choose >= 2*fmax/n/A
% fs_analog = fs_sym: analog sampling frequency (Hz), must satisfy fs_analog > 2*fmax
% n: length of chip sequences
% snr: signal-to-noise ratio (in dB)
% x: input string

% Output parameters:
% x_rec: sequence fed through DAD process



% This is the converter provided by Copilot with Matlab.
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

% Add white noise with SNR level snr.
% For reproducible noise samples (specify RNG seed)
% rng(0);                              % set seed
analog = analog + awgn(analog, snr, 'measured', 'db');

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

end

function omegaplot(plotset,snr,n)

% This function does log plots vs omegaw.

% Called by: main
% Calls: none

% Input parameters:
% n: length of chip sequences
% plotset: log(F) values (1st column cutoff 512, second column cutoff 1024)
% snr: signal-to-noise ratio (in dB)

% Plot using omegaw on x-axis
figure;
hold on;
x = (1:length(plotset))'; % x-axis (word rate)

h1 = plot(x, plotset(:,1), '-o', 'Color', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);
h2 = plot(x, plotset(:,2), '-^', 'Color', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

% Fit linear model y = m*x + b to both columns (since y already log-transformed)
ft = fittype('m*x + b', 'independent', 'x', 'coefficients', {'m','b'});
opts = fitoptions(ft);

colorsFit = [0 0.4470 0.7410; 0.8500 0.3250 0.0980]; % blue, orange
hfit = gobjects(2,1);
eqTexts = cell(2,1);

for col = 1:2
    y = plotset(:,col);
    % Start guesses: slope ~ (change in y)/(change in x), intercept ~ mean(y)
    dy = y(end) - y(1);
    dx = x(end) - x(1);
    m0 = dy / dx;
    b0 = mean(y) - m0 * mean(x);
    opts.StartPoint = [m0, b0];

    [cfun, gof] = fit(x, y, ft, opts);
    % Evaluate fit for smooth plotting
    xfit = linspace(min(x), max(x), 300)';
    yfit = cfun.m * xfit + cfun.b;
    hfit(col) = plot(xfit, yfit, '--', 'Color', colorsFit(col,:), 'LineWidth', 1.5);

    % Build LaTeX equation and R^2 string
    m = cfun.m;
    b = cfun.b;
    R2 = gof.rsquare;
    eqTexts{col} = sprintf('$\\log|F|(\\omega_{\\rm w}) = %.3g\\,\\omega_{\\rm w} %+.3g,\\; R^2 = %.4f$', m, b, R2);
end

% Place the two equations on the plot (upper-left and upper-right)
ax = gca;
xr = ax.XLim; yr = ax.YLim;
xpos2 = xr(1) + 0.45*(xr(2)-xr(1));
ypos2 = yr(2) +0.7*(yr(2)-yr(1));
text(30,3, eqTexts{1}, 'Interpreter', 'latex', 'FontSize', 12, 'BackgroundColor', 'none', 'Color', colorsFit(1,:));
text(25.5,5.2, eqTexts{2}, 'Interpreter', 'latex', 'FontSize', 12, 'BackgroundColor', 'none', 'Color', colorsFit(2,:));

% Create a string which represents the correct argument for snr:
if snr == Inf
    fstring = "$|F(\infty)|$";
else
    fstring = sprintf('$|F(%.3g)|$',snr);
end

% Then use the string in the title and y-axis.
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
tstring = append('log',fstring,' {\it vs}. word rate, ');
tstring = append(tstring,sprintf('$n=%d$',n));
title(tstring, 'Interpreter', 'latex');
ylabel(fstring);

% Labels, limits, legend
xlabel('Word rate $\omega_{\rm w}$ (Hz)');
xlim([1 length(plotset)]);
ylim([0 log(n)]);
grid on;
legend([h1 h2 hfit(1) hfit(2)], {'$\omega_a = 512$ Hz', '$\omega_a = 1024$ Hz', 'Exp fit (512)', 'Exp fit (1024)'}, ...
    'Location', 'best');

hold off;

end

function nomegaplot(plotset,snr,n,maxrun)

% This function does log plots vs omegaw in the case of noise.

% Called by: main
% Calls: none

% Input parameters:
% maxrun: number of simulated runs
% n: length of chip sequences
% plotset: log(F) values: 1st column one realization, then average, then intersection
% snr: signal-to-noise ratio (in dB)

% Plot using omegaw on x-axis
figure;
hold on;
x = (1:length(plotset))'; % x-axis (word rate)

h1 = plot(x, plotset(:,1), '-', 'Color', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);
h2 = plot(x, plotset(:,2), '-', 'Color', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);
h3 = plot(x, plotset(:,3), '-', 'Color', 'b', 'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

% Create a string which represents the correct argument for snr:
if snr == Inf
    fstring = "$|F(\infty)|$";
else
    fstring = sprintf('$|F(%.3g)|$',snr);
end

% Then use the string in the title and y-axis.
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
tstring = append('log',fstring,' {\it vs}. word rate, ');
tstring = append('log',fstring,sprintf(' {\\it vs}. word rate, $n=%d$, ',n));
tstring = append(tstring,sprintf('$\\omega_{\\rm a}=%d$ Hz, %d runs',n,maxrun));
title(tstring, 'Interpreter', 'latex');
ylabel(fstring);

% Labels, limits, legend
xlabel('Word rate $\omega_{\rm w}$ (Hz)');
xlim([1 length(plotset)]);
ylim([3 log(n)]);
grid on;
legend([h1 h2 h3], {'One Run', 'Average', 'Intersection'});

hold off;

end

function bandplot(plotset,snr,n,omegaw)

% This function does linear plots vs bandwidth.

% Called by: main
% Calls: none

% Input parameters:
% n: length of chip sequences
% plotset: 1st column cutoff omegaw, second column F
% snr: signal-to-noise ratio (in dB)

% Plot using omegaw on x-axis

% Then plot the result:
figure;
hold on;

% Column 1: black circles connected
h1 = plot(plotset(:,1),plotset(:,2), '-o', 'Color', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

% Fit linear model y = m*x + b to |F| vs bandwidth
x2 = plotset(:,1);
y2 = plotset(:,2);

ft2 = fittype('m*x + b', 'independent', 'x', 'coefficients', {'m','b'});
opts2 = fitoptions(ft2);
% Start guesses: slope ~ (dy/dx), intercept ~ mean(y)
dy2 = y2(end) - y2(1);
dx2 = x2(end) - x2(1);
m02 = dy2 / dx2;
b02 = mean(y2) - m02 * mean(x2);
opts2.StartPoint = [m02, b02];

[cfun2, gof2] = fit(x2, y2, ft2, opts2);

% Plot the linear fit line
xfit2 = linspace(min(x2), max(x2), 300)';
yfit2 = cfun2.m * xfit2 + cfun2.b;
hfit2 = plot(xfit2, yfit2, '--', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5);

% Build LaTeX equation string and display it on the plot
m = cfun2.m; b = cfun2.b; R2 = gof2.rsquare;
eqStr = sprintf('$|F|(\\omega_{\\rm a}) = %.3g\\,\\omega_{\\rm a} %+.3g,\\; R^2 = %.4f$', m, b, R2);

% Choose a location within axes to place the text (upper-left)
ax = gca;
xr = ax.XLim; yr = ax.YLim;
xpos = xr(1) + 0.05*(xr(2)-xr(1));
ypos = yr(2) - 0.08*(yr(2)-yr(1));
text(xpos, ypos, eqStr, 'Interpreter', 'latex', 'FontSize', 12, 'BackgroundColor', 'none', 'Color', [0.2 0.6 0.2]);

if snr == Inf
    fstring = "$|F(\infty)|$";
else
    fstring = sprintf('$|F(%.3g)|$',snr);
end

% Then use the string in the title and y-axis.
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
tstring = append(fstring,sprintf(' vs. bandwidth, $\\omega_{\\rm w}=%d$ Hz, ',omegaw));
tstring = append(tstring,sprintf('$n=%d$',n));
title(tstring, 'Interpreter', 'latex');
ylabel(fstring);

xlabel('Bandwidth $\omega_{\rm a}$ (Hz)');
xlim([512 1024]);
% ylim([0 n]);

% % Create legend labels "Analog cutoff = 512*column number" (use LaTeX)
% labels = arrayfun(@(c) sprintf('Analog cutoff $\omega_{\rm a}$ = %d$ Hz', 512*c), 1:2, 'UniformOutput', false);
% legend([h1 h2], labels, 'Location', 'best');
grid on;

hold off;

end

function nbandplot(plotset,snr,n,maxrun,omegaw)

% This function does regular plots vs bandwidth.

% Called by: main
% Calls: none

% Input parameters:
% n: length of chip sequences
% plotset: 1st column one realization, then average, then intersection
% snr: signal-to-noise ratio (in dB)

% Plot using omegaw on x-axis

% Then plot the result:
figure;
hold on;

h1 = plot(plotset(:,4), plotset(:,1), '-', 'Color', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);
h2 = plot(plotset(:,4), plotset(:,2), '-', 'Color', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);
h3 = plot(plotset(:,4), plotset(:,3), '-', 'Color', 'b', 'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

% Create a string which represents the correct argument for snr:
if snr == Inf
    fstring = "$|F(\infty)|$";
else
    fstring = sprintf('$|F(%.3g)|$',snr);
end

% Then use the string in the title and y-axis:
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
tstring = append(fstring,sprintf(' vs. bandwidth, $\\omega_{\\rm w}=%d$ Hz, ',omegaw));
tstring = append(tstring,sprintf('$n=%d$, %d runs',n,maxrun));
title(tstring, 'Interpreter', 'latex');
ylabel(fstring);

xlabel('Bandwidth $\omega_{\rm a}$ (Hz)');
xlim([512 1024]);
% ylim([0 n]);

% % Create legend labels "Analog cutoff = 512*column number" (use LaTeX)
% labels = arrayfun(@(c) sprintf('Analog cutoff $\omega_{\rm a}$ = %d$ Hz', 512*c), 1:2, 'UniformOutput', false);
% legend([h1 h2], labels, 'Location', 'best');
grid on;
legend([h1 h2 h3], {'One Run', 'Average', 'Intersection'},'Location', 'northwest');

hold off;

end

function faithful = getfaithful(n,omegaw,fmax,snr)
% This function calculates the faithful set of Walsh matrix rows given a
% signal-noise ratio.

% Called by: main
% Calls: bestmatch

% Input parameters:
% fmax: desired analog bandwidth (Hz)
% n: size of Walsh matrix
% omegaw: word transmission rate.
% snr: signal-to-noise ratio (dB)

% Output variables:
% bestmatch: index of encoding which best matches the received state

% Internal variables:
% DADflag: 1 if using Matlab's converter; 2 if using Adam's

% Compute the Walsh matrix of size n.
h = n*fwht(eye(n));

for i = 1:n
    % See if this can be vectorized!
    % Compute the best match for row i.
    bestvec(i) = bestmatch(i,h,n,omegaw,fmax,snr);
end

% Echo value as the algorithm progresses.
if mod(omegaw,8)==0
    fprintf('omegaw = %d\n', omegaw);
end

% Compute the faithful set (where bestvec (output) = input row):
faithful = (bestvec == (1:length(bestvec)));

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

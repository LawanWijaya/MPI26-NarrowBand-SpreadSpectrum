%{
daefadeplot.m

This code performs several algorithms to take the best codebook from
the feasible set.

This code uses the Signal Processing Toolbox, 
Communication Toolbox and the Curve Fitting Toolbox.

Last modified by David A. Edwards on 6/17/26.

%}

% Clear all variables from previous runs.
clear all


% Define parameters.

fmax = 512; % desired analog bandwidth (Hz): also called omegaa.
% i: index of rows
n = 512; % Size of Hadamard matrix
% omegaw: Word transmission rate.
omegawmax = 64; % Maximum word transmission rate.
w = 4; % Word length

% Variables
bestvec = zeros(n,1); % Vector of row of best matches.
% faithful: Vector of faithful set.
plotset = zeros(omegawmax,2); % Matrix of plotting vectors.

% % Set LaTeX interpreters and default font size for all plots.
set(groot, 'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesFontSize', 14, ...
    'defaultTextFontSize', 14);

omegaw = 10;

% Construct the faithful matrices.

% Number or runs to use with noise.
maxrun = 50;
maxrun = 10;
% rng(0);

faderun = 100;

pc = 0.9449;  % probability of staying correct
pe  = 0.8509;  % probability of propagating error

% Define plot arrays.
% faithmat contains the results from all the runs.  We use a cell array so
% that we can store the arrays inside without worrying about three indices:
faithmat = cell(maxrun,1);
plotset = zeros(omegawmax,3);

% Compute the Walsh matrix of size n.
h = n*fwht(eye(n));
eset = 2^w;

snr = -3;

% for snr = -3:-6:-9

for run = 1:maxrun
    faithful = getfaithful(h,omegaw,fmax,snr);
    faithmat{run,1} = faithful;
end

% Adam wants the greedy algorithm alone given fcap, so:

% Compute the feasible set.
% Step 1: Combine the rows of faithmat into a single matrix which can be
% analyzed.
% faithmat is maxrun-by-1 cell, each cell is a 1-by-N or N-by-1 vector of 0/1
onemat = vertcat(faithmat{:});        % maxrun-by-N numeric matrix

% Step 2: Create an array that has a 1 where all the rows where EVERY entry
% is 1:
twomat = all(onemat==1, 1)'; % 1-by-N logical: true where every run has a 1

% Step 3: Now create an array that lists the row number (Walsh score) of
% every matching row.
feasible = find(twomat);

greedy = codeset(feasible,w);
% Generate the codebooks with each algorithm.
codebook = {h(greedy(:,1),:),h(greedy(:,2),:)};



% % First plot: Heat map to show which code words transmit best.
%
% heatplot = cell(1,2);
%
% % Have to check each algorithm.
% for alg = 1:2
%     heatplot{alg}=zeros(eset,eset);
%     for frun = 1:faderun
%         % Need to reset this each time.
%     tempmat = zeros(eset,eset);
%     % Get the received signal.  Note that we are passing the entire set
%     % of rows, so
%     faded = simulate_fading(codebook{alg},pc,pe);
%     % Then compute the Walsh score of the input and output, which is 1 more than the number of
%     % changes.
%     outwalsh = sum(faded(:,1:end-1) .* faded(:,2:end) == -1, 2 );
%     % [greedy(:,alg),outwalsh]
%     % Construct a matrix of the distances between greedy(i) and outwalsh(i):
%     D = abs(greedy(:,alg) - outwalsh.');
%     % Then find the minimum entry in each row, and find the index that matches
%     % it.
%     [distPerRow, cols] = min(D, [], 1);
%     % Then if cols is the same as the input row number, we've matched.
%     ind = sub2ind([eset,eset], 1:eset, cols(:).');      % row and col subscripts are both 1-by-n
%     tempmat(ind) = 1;
%     % Now tempmat has 1s only where they match, so you add it.
%     heatplot{alg} = heatplot{alg} + tempmat;
%     % [1:16;cols]'
%     end
%     heatplot{alg} = heatplot{alg}/faderun;
%     corrplot(heatplot{alg},snr,n,maxrun,omegaw,faderun,alg,pc,pe)
% end
%
% % Second plot: Heat map to show how the system is affected when the
% % probabilities change.  Note that this only affects the fading, so we
% % can use the same codebook.
%
% % Now do the computation for the probability map.
% probmat = linspace(0.8,1,20);
%
% hp2 = cell(1,2);
%
% for alg = 1:2
%     hp2{alg}=zeros(20,20);
%     for i = 1:length(probmat)
%         for j = 1:length(probmat)
%             for frun = 1:faderun
%                 % Get the received signal.  Note that we are passing the entire set
%                 % of rows, so
%                 faded = simulate_fading(codebook{alg},probmat(i),probmat(j));
%                 % Then compute the Walsh score of the input and output, which is 1 more than the number of
%                 % changes.
%                 outwalsh = sum(faded(:,1:end-1) .* faded(:,2:end) == -1, 2 );
%                 % [greedy(:,alg),outwalsh]
%                 % Construct a matrix of the distances between greedy(i) and outwalsh(i):
%                 D = abs(greedy(:,alg) - outwalsh.');
%                 % Then find the minimum entry in each row, and find the index that matches
%                 % it.
%                 [distPerRow, cols] = min(D, [], 1);
%                 % % Then if cols is the same as the input row number, we've matched.
%                 % So we just sum them up.
%                 hp2{alg}(i,j)=hp2{alg}(i,j)+sum((1:16)==cols);
%                 % [1:16;cols;(1:16)==cols]'
%             end
%         end
%     end
%     % Then we normalize by the number of runs and the total length:
%     hp2{alg} = hp2{alg}/faderun/eset;
%     corrbyp(hp2{alg},snr,n,maxrun,omegaw,faderun,alg,probmat)
% end

% % Third plot: probability vs bandwidth, with noise.
% % for snr = -3:-6:-9
% for snr = -3:-3
%     % Reset the row number.
%     count = 1;
%     for fmax = 512:8:1024 % Go from 512 to 1024 by 8s
%         for run = 1:maxrun
%             faithful = getfaithful(h,omegaw,fmax,snr);
%             faithmat{run,1} = faithful;
%         end
% 
%         % Echo value
%         if mod(fmax,64)==0
%             fprintf('fmax = %d\n', fmax);
%         end
% 
%         % Adam wants the greedy algorithm alone given fcap, so:
% 
%         % Compute the feasible set.
%         % Step 1: Combine the rows of faithmat into a single matrix which can be
%         % analyzed.
%         % faithmat is maxrun-by-1 cell, each cell is a 1-by-N or N-by-1 vector of 0/1
%         onemat = vertcat(faithmat{:});        % maxrun-by-N numeric matrix
% 
%         % Step 2: Create an array that has a 1 where all the rows where EVERY entry
%         % is 1:
%         twomat = all(onemat==1, 1)'; % 1-by-N logical: true where every run has a 1
% 
%         % Step 3: Now create an array that lists the row number (Walsh score) of
%         % every matching row.
%         feasible = find(twomat);
% 
%         greedy = codeset(feasible,w);
%         % Generate the codebooks with each algorithm.
%         codebook = {h(greedy(:,1),:),h(greedy(:,2),:)};
% 
%         for alg = 1:2
%             hitrate = 0;
%             for frun = 1:faderun
%                 % Get the received signal.  Note that we are passing the entire set
%                 % of rows, so
%                 faded = simulate_fading(codebook{alg},0.99,0.99);
%                 % Then compute the Walsh score of the input and output, which is 1 more than the number of
%                 % changes.
%                 outwalsh = sum(faded(:,1:end-1) .* faded(:,2:end) == -1, 2 );
%                 % [greedy(:,alg),outwalsh]
%                 % Construct a matrix of the distances between greedy(i) and outwalsh(i):
%                 D = abs(greedy(:,alg) - outwalsh.');
%                 % Then find the minimum entry in each row, and find the index that matches
%                 % it.
%                 [distPerRow, cols] = min(D, [], 1);
%                 % % Then if cols is the same as the input row number, we've matched.
%                 % So we just sum them up.
%                 hitrate = hitrate +sum((1:16)==cols);
%                 % [1:16;cols;(1:16)==cols]'
%             end
%             plotset(count,alg+1) = hitrate/faderun/eset;
%         end
%         % We put the x-value in the first column for plotting.
%         plotset(count,1) = fmax;
%         count = count+1;
% 
%     end
% 
% end
% 
% % Then plot the result:
% 
% bandplot(plotset,snr,n,maxrun,omegaw,faderun,0.99,0.99);

% Fourth plot: loop over omegaw, add noise, but just for n=512.

% Loop over snr values.

% for snr = -3:-6:-9

plotset = zeros(omegawmax,3);

for snr = -3:-3
    for omegaw = 1:omegawmax
        if mod(omegaw,8)==0
            fprintf('omegaw = %d\n', omegaw);
        end

        for run = 1:maxrun
            faithful = getfaithful(h,omegaw,fmax,snr);
            faithmat{run,1} = faithful;
        end

        % Adam wants the greedy algorithm alone given fcap, so:

        % Compute the feasible set.
        % Step 1: Combine the rows of faithmat into a single matrix which can be
        % analyzed.
        % faithmat is maxrun-by-1 cell, each cell is a 1-by-N or N-by-1 vector of 0/1
        onemat = vertcat(faithmat{:});        % maxrun-by-N numeric matrix

        % Step 2: Create an array that has a 1 where all the rows where EVERY entry
        % is 1:
        twomat = all(onemat==1, 1)'; % 1-by-N logical: true where every run has a 1

        % Step 3: Now create an array that lists the row number (Walsh score) of
        % every matching row.
        feasible = find(twomat);

        greedy = codeset(feasible,w);

        codebook = {h(greedy(:,1),:),h(greedy(:,2),:)};

        for alg = 1:2
            hitrate = 0;
            for frun = 1:faderun
                % Get the received signal.  Note that we are passing the entire set
                % of rows, so
                faded = simulate_fading(codebook{alg},0.99,0.99);
                % Then compute the Walsh score of the input and output, which is 1 more than the number of
                % changes.
                outwalsh = sum(faded(:,1:end-1) .* faded(:,2:end) == -1, 2 );
                % [greedy(:,alg),outwalsh]
                % Construct a matrix of the distances between greedy(i) and outwalsh(i):
                D = abs(greedy(:,alg) - outwalsh.');
                % Then find the minimum entry in each row, and find the index that matches
                % it.
                [distPerRow, cols] = min(D, [], 1);
                % % Then if cols is the same as the input row number, we've matched.
                % So we just sum them up.
                hitrate = hitrate +sum((1:16)==cols);
                % [1:16;cols;(1:16)==cols]'
            end
            plotset(omegaw,alg+1) = hitrate/faderun/eset;
        end
        % We put the x-value in the first column for plotting.
        plotset(:,1) = 1:omegawmax;
    end

    % % Then plot the result:

    omegaplot(plotset,snr,n,maxrun,faderun,0.99,0.99);

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

function corrplot(permat,snr,n,maxrun,omegaw,faderun,alg,pc,pe)

% This function does plots vs omegaw of the minimum distance.

% Called by: main
% Calls: none

% Input parameters:
% maxrun: number of simulated runs
% n: length of chip sequences
% plotset: log(F) values (1st column cutoff 512, second column cutoff 1024)
% snr: signal-to-noise ratio (in dB)

if alg==1
    line1 = 'Heat Map for Greedy First, Percent Successfully Decoded';
else
    line1 = 'Heat Map for Greedy Last, Percent Successfully Decoded';
end

% Plot using omegaw on x-axis
figure;
% % Then generate the heat map.  Here we must set the h or things won't
% work:
h = heatmap(permat, 'Colormap', parula);
h.ColorbarVisible = 'on';
h.CellLabelFormat = '%.2g';    % format the numbers shown in cells
h.YLabel = 'Transmitted';
h.XLabel = 'Decoded';
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
line1 = append(line1,sprintf(', $n=%d$, $\\omega_{\\rm a}=%d$ Hz', n,n));
line2 = sprintf('$\\omega_{\\rm w}=%d$ Hz, %d AWGN runs, %d fade runs, $r=%d$', omegaw, maxrun, faderun, snr);
line2 = append(line2,sprintf(', $P_{\\rm c}=%.4f$, $P_{\\rm e}=%.4f$',pc,pe));


% Put each line in a cell array so title creates multiple lines
h.Title = {line1; line2};
h.Interpreter = 'latex';

end

function omegaplot(plotset,snr,n,maxrun,faderun,pc,pe)

% This function does plots vs omegaw of the minimum distance.

% Called by: main
% Calls: none

% Input parameters:
% maxrun: number of simulated runs
% n: length of chip sequences
% plotset: log(F) values (1st column cutoff 512, second column cutoff 1024)
% snr: signal-to-noise ratio (in dB)

% Plot using omegaw on x-axis
figure;
hold on;
% x = (1:length(plotset))'; % x-axis (word rate)

h1 = plot(plotset(:,1), plotset(:,2), '-', 'Color', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);
h2 = plot(plotset(:,1), plotset(:,3), '-', 'Color', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

% Then use the string in the title and y-axis:
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
line1 = 'Probability correct {\it vs}. word rate';
line1 = append(line1,sprintf(', $n=%d$, $\\omega_{\\rm a}=%d$ Hz', n,n));
line2 = sprintf('%d AWGN runs, %d fade runs, $r=%d$', maxrun, faderun, snr);
line2 = append(line2,sprintf(', $P_{\\rm c}=%.4f$, $P_{\\rm e}=%.4f$',pc,pe));
tstring = {line1;line2};
title(tstring);
ylabel('Probability correct');

% Labels, limits, legend
xlabel('Word rate $\omega_{\rm w}$ (Hz)');
xlim([1 length(plotset)]);
grid on;
legend([h1 h2], {'Greedy First', 'Greedy Last'});

hold off;

end

function corrbyp(permat,snr,n,maxrun,omegaw,faderun,alg,probmat)

% This function does plots vs omegaw of the minimum distance.

% Called by: main
% Calls: none

% Input parameters:
% maxrun: number of simulated runs
% n: length of chip sequences
% plotset: log(F) values (1st column cutoff 512, second column cutoff 1024)
% snr: signal-to-noise ratio (in dB)

if alg==1
    line1 = 'Heat Map for Greedy First, Percent Successfully Decoded';
else
    line1 = 'Heat Map for Greedy Last, Percent Successfully Decoded';
end

% Plot using omegaw on x-axis
figure;
% % Then generate the heat map.  Here we must set the h or things won't
% work:
% h = heatmap(permat);
h = heatmap(permat, 'Colormap', parula);
% Always put this line first!
h.Interpreter = 'latex';
h.ColorbarVisible = 'on';
h.CellLabelFormat = '%.2g';    % format the numbers shown in cells
% DON'T USE THIS FORM!
% h.XDisplayData = probmat;             % x axis values (1-by-c or c-by-1)
% h.YDisplayData = probmat;             % y axis values (1-by-r or r-by-1)

% If you want custom labels (e.g. formatted strings):
h.XDisplayLabels = compose('%.2f', probmat);   % cell/string array of labels (c elements)
h.YDisplayLabels = compose('%.2f', probmat);   % cell/string array of labels (c elements)
h.XLabel = '$P_{\rm c}$';
h.YLabel = '$P_{\rm e}$';
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
line1 = append(line1,sprintf(', $n=%d$, $\\omega_{\\rm a}=%d$ Hz', n,n));
line2 = sprintf('$\\omega_{\\rm w}=%d$ Hz, %d AWGN runs, %d fade runs, $r=%d$', omegaw, maxrun, faderun, snr);
% line2 = append(line2,sprintf(', $P_{\\rm c}=%.4f$, $P_{\\rm e}=%.4f$',pc,pe));


% Put each line in a cell array so title creates multiple lines
h.Title = {line1; line2};

end

function bandplot(plotset,snr,n,maxrun,omegaw,faderun,pc,pe)

% This function does linear plots vs bandwidth.

% Called by: main
% Calls: none

% Input parameters:
% n: length of chip sequences
% plotset: plot of probabilities
% snr: signal-to-noise ratio (in dB)

% Plot using omegaw on x-axis

% Then plot the result:
figure;
hold on;

h1 = plot(plotset(:,1), plotset(:,2), '-', 'Color', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);
h2 = plot(plotset(:,1), plotset(:,3), '-', 'Color', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

% Then use the string in the title and y-axis:
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
line1 = 'Probability correct {\it vs}. bandwidth';
line1 = append(line1,sprintf(', $n=%d$, $\\omega_{\\rm w}=%d$ Hz', n,omegaw));
line2 = sprintf('%d AWGN runs, %d fade runs, $r=%d$', maxrun, faderun, snr);
line2 = append(line2,sprintf(', $P_{\\rm c}=%.4f$, $P_{\\rm e}=%.4f$',pc,pe));
tstring = {line1;line2};
title(tstring);
ylabel('Probability correct');

xlabel('Bandwidth $\omega_{\rm a}$ (Hz)');
xlim([512 1024]);
% ylim([0 n]);

grid on;
legend([h1 h2], {'Greedy First', 'Greedy Last'},'Location', 'northwest');

hold off;

end

function faithful = getfaithful(h,omegaw,fmax,snr)
% This function calculates the faithful set of Walsh matrix rows given a
% signal-noise ratio.  This version is optimized so h doesn't have to be
% computed each time.

% Called by: main
% Calls: bestmatch

% Input parameters:
% fmax: desired analog bandwidth (Hz)
% h: Walsh matrix
% omegaw: word transmission rate.
% snr: signal-to-noise ratio (dB)

% Output variables:
% bestmatch: index of encoding which best matches the received state

% Internal variables:
% DADflag: 1 if using Matlab's converter; 2 if using Adam's

n = size(h,2); % size of Walsh matrix

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

function greedy = codeset(feasible,w)

% This function computes the codeset using two different greedy algorithms.

% Called by: main
% Calls: none

% Input variables:
% feasible: vector of indices of faithful codes
% w: length of word

% Output variable:
% greedy: two columns of codesets generated by different greedy algorithms

% Internal variables:
esize = 2^w; % size of codebook

greedy = zeros(esize,2); % Matrix of encoding set

% feasible(end) is the largest Walsh score, so the optimal distance would
% be this over the size of the encoding set, no matter the algorithm
greedy(esize,:)=feasible(end);

% Now do the greedy first algorithm.
for j = esize:-1:2
    % Step 1.  Set the optimum distance for what is left, using the basin
    % of attraction distance.
    optdis = greedy(j,1)/(j-1/2);
    % Step 2.  Filter the feasible vectors so only entries are retained
    % which are an odd distance from the set vectors.  First, construct a
    % vector of 0's and 1's depending on whether the distance is odd:
    indset = mod(abs(feasible - greedy(j,1)),2)==1;
    % Then select only the entries of feasible that are spaced correctly:
    oddfeas = feasible(indset);
    % Step 3.  Choose the entry that is optdis away or more.  Find finds a
    % vector of all the INDICES that satsify, and then we choose the last
    % (highest) one.
    testfind = find(oddfeas <= (greedy(j,1)-optdis),1,'last');
    testfind = oddfeas(testfind);
    % If testfind can't find anything (probably only at the final
    % point), then we just reduce the optimal distance by 1 until we
    % get a result.
    while isempty(testfind)
        optdis = optdis - 1;
        testfind = find(oddfeas <= (greedy(j,1)-optdis),1,'last');
        testfind = oddfeas(testfind);
    end
    greedy(j-1,1) = testfind;
end

% Now do the greedy last algorithm.
for j = esize:-1:2
    % Step 1.  Set the optimum distance for what is left, using the basin
    % of attraction distance.
    optdis = greedy(j,2)/(j-1/2);
    % Step 2.  Filter the feasible vectors so only entries are retained
    % which are an odd distance from the set vectors.  First, construct a
    % vector of 0's and 1's depending on whether the distance is odd:
    indset = mod(abs(feasible - greedy(j,2)),2)==1;
    % Then select only the entries of feasible that are spaced correctly:
    oddfeas = feasible(indset);
    % Step 3.  Choose the entry that is optdis away or more.  Find finds a
    % vector of all the INDICES that satsify, and then we choose the last
    % (highest) one.
    testfind = find(oddfeas <= (greedy(j,2)-optdis),1,'last');
    testfind = oddfeas(testfind);
    % If testfind can't find anything (probably only at the final
    % point), then we just reduce the optimal distance by 1 until we
    % get a result.
    while isempty(testfind)
        optdis = optdis - 1;
        testfind = find(oddfeas >= (greedy(j,2)-optdis),1,'first');
        testfind = oddfeas(testfind);
    end
    greedy(j-1,2) = testfind;
end

end

function return_data = simulate_fading(inpseq,pc,pe)
% SIMULATE_FADING functions call take two forms: either submit only an
% input sequence, and default probabilities will be applied, or manually
% specify the probabilities.  To do that, you do the function call
% simulated_fading(inpseq,"stay_good",prob,"stay_bad",prob);
arguments (Input)
    inpseq                  % the input sequence
    pc % probability of staying correct
    pe % probability of propagating error
end
% This block says that the values for the probabilities are optional.  If I
% don't include them, the code will fill in these values.
% arguments (Input)
%     options.stay_good = 0.9449  % probability of staying faithful
%     options.stay_bad  = 0.8509  % probability of propagating error
% end

% Internal variables:
% i: row number in codebook
% j: chip number in code

% State variable.  If zero, then the bit sequence is correct.  If 1, it has
% an error.
state = 0;
outseq = inpseq;

% rng(0);

[imax,jmax]=size(inpseq);

for i = 1:imax % Index over rows
    for j = 1:length(inpseq) % Index over chips
        if state==1 % preceding i has error
            if rand > pe   % error ends
                state = 0;
                % disp(j)
            else                         % error persists
                outseq(i,j) = -1*inpseq(i,j);
            end
        else     % preceding i is correct
            if rand > pc  % begin error
                % disp(j);
                state = 1;
                outseq(i,j) = -1*inpseq(i,j);
            end
        end
    end
end

return_data = outseq;

end
%{
daeintraplot.m

This code performs several algorithms to take the best codebook from
the feasible set and generates plots.

This code uses the Signal Processing Toolbox, 
Communication Toolbox and the Curve Fitting Toolbox.

Last modified by David A. Edwards on 6/17/26.

%}


% Define parameters.

fmax = 512; % desired analog bandwidth (Hz): also called omegaa.
% i: index of rows
n = 512; % Size of Hadamard matrix
% omegaw: Word transmission rate.
omegawmax = 64; % Maximum word transmission rate.
w = 4; % word length

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

% Construct the faithful matrices.

% Number or runs to use with noise.
maxrun = 50;
maxrun = 10;

% Define plot arrays.
% faithmat contains the results from all the runs.  We use a cell array so
% that we can store the arrays inside without worrying about three indices:
faithmat = cell(maxrun,1);
plotset = zeros(omegawmax,3);

% First plot: loop over omegaw, add noise, but just for n=512.

% Loop over snr values.

for snr = -3:-6:-9

    for omegaw = 1:omegawmax

        for run = 1:maxrun
            faithful = getfaithful(n,omegaw,fmax,snr);
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

        % Then get the minimum distance in each column for plotting.  diff gives us
        % the differences, and the other arguments in min tell it to minimize by
        % columns.
        % We put the x-value in the first column for plotting.
        plotset(omegaw,:) = [omegaw,min(diff(greedy,1,1), [], 1)];  % 1-by-numCols

    end

    % % Then plot the result:

    omegaplot(plotset,snr,n,maxrun);

end

omegaw = 10;
% Second plot: F vs bandwidth, with noise.
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

        % Then get the minimum distance in each column for plotting.  diff gives us
        % the differences, and the other arguments in min tell it to minimize by
        % columns.
        % We put the x-value in the first column for plotting.
        plotset(count,:) = [fmax,min(diff(greedy,1,1), [], 1)];  % 1-by-numCols
        count = count+1;

    end

    % % Then plot the result:

    bandplot(plotset,snr,n,maxrun,omegaw);

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

function omegaplot(plotset,snr,n,maxrun,omegaw)

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

% % Create a string which represents the correct argument for snr:
% if snr == Inf
%     fstring = "$|F(\infty)|$";
% else
%     fstring = sprintf('$|F(%.3g)|$',snr);
% end

% Then use the string in the title and y-axis:
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
tstring = sprintf('Minimum distance {\\it vs}. word rate, $n=%d$, ',n);
tstring = append(tstring,sprintf('$\\omega_{\\rm a}=%d$ Hz, %d runs, $r=%d$',n,maxrun,snr));
title(tstring, 'Interpreter', 'latex');
ylabel('Minimum distance');

% Labels, limits, legend
xlabel('Word rate $\omega_{\rm w}$ (Hz)');
xlim([1 length(plotset)]);
grid on;
legend([h1 h2], {'Greedy First', 'Greedy Last'});

hold off;

end

function bandplot(plotset,snr,n,maxrun,omegaw)

% This function does linear plots vs bandwidth.

% Called by: main
% Calls: none

% Input parameters:
% n: length of chip sequences
% plotset: log(F) values (1st column cutoff 512, second column cutoff 1024)
% snr: signal-to-noise ratio (in dB)

% Plot using omegaw on x-axis

% Then plot the result:
figure;
hold on;

h1 = plot(plotset(:,1), plotset(:,2), '-', 'Color', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);
h2 = plot(plotset(:,1), plotset(:,3), '-', 'Color', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 6);

% % Create a string which represents the correct argument for snr:
% if snr == Inf
%     fstring = "$|F(\infty)|$";
% else
%     fstring = sprintf('$|F(%.3g)|$',snr);
% end

% Then use the string in the title and y-axis:
% IMPORTANT: If you have a STRING, use typical LaTEX notation.
%               If you have SPRINTF, then use \\ everywhere.
tstring = 'Minimum distance {\it vs}. bandwidth, ';
tstring = append(tstring,sprintf('$\\omega_{\\rm w}=%d$ Hz, $n=%d$, ',omegaw,n));
tstring = append(tstring,sprintf('%d runs, $r=%d$',maxrun,snr));
title(tstring);
ylabel('Minimum distance');

xlabel('Bandwidth $\omega_{\rm a}$ (Hz)');
xlim([512 1024]);
% ylim([0 n]);

grid on;
legend([h1 h2], {'Greedy First', 'Greedy Last'},'Location', 'northwest');

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

% Sort the columns in increasing order.
greedy = sort(greedy);

end
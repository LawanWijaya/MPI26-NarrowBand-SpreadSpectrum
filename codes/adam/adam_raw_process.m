%% Boilerplate

% Add all functions in project to path
addpath(genpath('C:\Users\adamn\OneDrive\Documents\GitHub\MPI26-NarrowBand-SpreadSpectrum'))

%% Compare algorithms for constructing codebooks
% The primary deliverable of this subblock is a figure

% Set parameters
wr = 4;   % word size (ie bits/word)
n = 9;    % code size (ie chips/code)
cf = 700; % cutoff for low-pass filter (interpreted in Hz)

% Collect set of faithful rows to test with
[faithful_coarse,faithful_fine] = run_cutoff_trend(n,cf);
faithful_set = find(faithful_fine);

% Collection algorithms to test
algos = {@lawan5_best @david1_best @david2_best @henry1_best};
algo_titles = ["Lawan" "David1" "David2" "Henry"];

% Instantiate figure for visualization
algo_comparison = figure(6);
clf(algo_comparison);
tiles = tiledlayout(2,2,'Padding','loose','TileSpacing','compact');
ax = gobjects(4,1); % WARNING: the numbers here (1,2,4) are all chosen to
                    % to fit the total number of algorithms. To test a
                    % longer (or shorter) list of algorithms, one would
                    % need to adjust these values

% Maximum basin size (for standardizing separate plots in figure)
maxbs = 0;

for k = 1:4

    % Move to next plot in figure
    ax(k) = nexttile;

    % Prepare for multiple graphics
    hold on

    % Construct collection of codebooks
    lawbook = construct_lawbook(wr,faithful_set,algos{k});
    cd_tot = size(lawbook,1);

    % Compute minimum basin sizes for each codebook in a collection
    minbs = zeros(cd_tot,1);
    for i = 1:cd_tot
        minbs(i) = test_basin(lawbook(i,:));
    end

    % Update maximum basin size
    maxbs = max(maxbs,max(minbs));

    % Prepare visual representation of codeboks by constructing
    % 'checkerboard' pattern: 0 if row is included, 1 otherwise
    code_swapping = zeros(cd_tot,max(faithful_set));
    for i = 1:cd_tot
        code_swapping(i,lawbook(i,:)) = true;
    end
    code_swapping = 1-code_swapping;

    % Represent codebook entries for current algorithm
    imagesc(ax(k),1:cd_tot,1:max(faithful_set),code_swapping.');
    colormap([0 0 0; 1 1 1]);

    % Present minimum basin sizes for codebooks constructed with current
    % algorithm
    x = (1:cd_tot);
    yyaxis(ax(k),'right')
    plot(ax(k),x,minbs,'r-o','LineWidth',2)

    % Title current algorithm
    title(ax(k),sprintf("Codebooks by %s",algo_titles(k)))

end

% Format axes for consistency
for k = 1:4
    yyaxis(ax(k),'right')
    ax(k).YColor = 'k';
    ylim(ax(k),[0,maxbs+1])
    xlim(ax(k),[0,cd_tot+1])
end

% Place figure
tiles.OuterPosition = [0.05 0.05 0.85 0.9];

% Label axes
linkaxes(ax,'x')
title(tiles,"Generation of Codebooks")
xlabel(tiles,"Codebook #")
ylabel(tiles,"Code (Walsh) #")
annotation('textbox',[0.88 0.57 0.1 0.1],...
    'String','Min Basin Size',...
    'Rotation',-90,...
    'EdgeColor','none',...
    'HorizontalAlignment','center', ...
    'FontSize',12);

% Save figure
saveas(tiles,'comparing_algos.png')

%% Visualize faith as function of cutoff
% The primary deliverable of this subblock is a figure

% Set parameters
cutoffs = 512:2:1024;   % range of cutoffs to test
n = 9;                  % number of chips in an encoding

% Determine faithful sets
[counts,data] = run_cutoff_trend(n,cutoffs);

% Set up figure and plot 'checkerboard': 1 if row is faithful, 0 otherwise
checkerboard= figure(3);
clf(checkerboard);
checkerboard_ax = axes(checkerboard);
imagesc(checkerboard_ax,cutoffs,1:(2^n),data);
colormap([0 0 0; 1 1 1]);

% Label axes
xlabel(checkerboard_ax,'Cutoff (Bandwidth) in Hz')
ylabel(checkerboard_ax,'Walsh Row #')
title(checkerboard_ax,'Trends in Faithfulness')

% Format colorbar
cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = {'Unfaithful','Faithful'};
cb.Position = [0.75 0.7 0.02 0.15];

% Save figure
saveas(checkerboard_ax,'checkerboard.png')
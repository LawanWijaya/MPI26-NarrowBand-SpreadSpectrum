
addpath(genpath('C:\Users\adamn\OneDrive\Documents\GitHub\MPI26-NarrowBand-SpreadSpectrum'))

%%

[faithful_coarse,faithful_fine] = run_cutoff_trend(9,700);
faithful_set = find(faithful_fine);

length(faithful_set)

wr = 4;
n = 9;

lawbook = construct_lawbook(wr,faithful_set,@david1_best);
lawbook(1,:);

basins = zeros(size(lawbook,1),1);
for i = 1%:size(lawbook,1)
    basins(i) = test_basin(lawbook(i,:));
end

lawbook(1,:);
%basins

basins;

code_swapping = zeros(size(lawbook,1),max(faithful_set));

for i = 1:size(lawbook,1)
    code_swapping(i,lawbook(i,:)) = true;
end
code_swapping = 1-code_swapping;

checkerboard_s = figure(4);
clf(checkerboard_s);
checkerboard_s_ax = axes(checkerboard_s);
imagesc(checkerboard_s_ax,1:size(lawbook,1),1:(2^n),code_swapping.');
colormap([0 0 0; 1 1 1]);

xlabel(checkerboard_s_ax,'Codebok #')
ylabel(checkerboard_s_ax,'Walsh Row #')
title(checkerboard_s_ax,'Lawbook')

saveas(checkerboard_s_ax,'checkerboard_s.png')

%%

cutoffs = 512:2:1024;

n = 9;

[counts,data] = run_cutoff_trend(n,cutoffs);

%%

% visualize cont output
trend = figure(1);
clf(trend);
trend_ax = axes(trend);
plot(trend_ax,cutoffs,counts,'o');
title(trend_ax,'Faithful Codes');

%%

checkerboard= figure(3);
clf(checkerboard);
checkerboard_ax = axes(checkerboard);
imagesc(checkerboard_ax,cutoffs,1:(2^n),data);
colormap([0 0 0; 1 1 1]);

xlabel(checkerboard_ax,'Cutoff (Bandwidth) in Hz')
ylabel(checkerboard_ax,'Walsh Row #')
title(checkerboard_ax,'Trends in Faithfulness')

cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = {'Unfaithful','Faithful'};
cb.Position = [0.75 0.7 0.02 0.15];

saveas(checkerboard_ax,'checkerboard.png')
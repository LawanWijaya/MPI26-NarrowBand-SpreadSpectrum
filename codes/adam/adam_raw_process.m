
addpath(genpath('C:\Users\adamn\OneDrive\Documents\GitHub\MPI26-NarrowBand-SpreadSpectrum'))

%%
cutoffs = 600:30:1500;

[counts,data] = run_cutoff_trend(512,cutoffs);

%%

% visualize cont output
trend = figure(1);
clf(trend);
trend_ax = axes(trend);
plot(trend_ax,cutoffs,counts,'o');
title(trend_ax,'Faithful Codes');

%%

checkerboard= figure(2);
clf(checkerboard);
checkerboard_ax = axes(checkerboard);
imagesc(checkerboard_ax,data);
colormap([0 0 0; 1 1 1]);
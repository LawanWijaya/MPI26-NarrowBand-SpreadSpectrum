
cutoffs = 600:30:1500;

counts = run_cutoff_trend(512,cutoffs);

% visualize cont output
trend = figure(2);
clf(trend);
trend_ax = axes(trend);
plot(trend_ax,cutoffs,counts,'o');
title(trend_ax,'Faithful Codes');
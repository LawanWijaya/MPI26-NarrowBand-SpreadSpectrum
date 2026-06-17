
addpath(genpath('C:\Users\adamn\OneDrive\Documents\GitHub\MPI26-NarrowBand-SpreadSpectrum'))

%%

[faithful_coarse,faithful_fine] = run_cutoff_trend(9,700);
faithful_set = find(faithful_fine);

wr = 4;
n = 9;

algos = {@lawan5_best @david1_best @david2_best @henry1_best};
algo_titles = ["Lawan" "David" "Dr. Edwards" "Henry"];

algo_comparison = figure(6);
clf(algo_comparison);

t = tiledlayout(2,2,'Padding','loose','TileSpacing','compact');
ax = gobjects(4,1);

M_temp = 0;

for k = 1:4

    ax(k) = nexttile;

    lawbook = construct_lawbook(wr,faithful_set,algos{k});
    cd_tot = size(lawbook,1);

    lawbook(1,end)

    basins = zeros(cd_tot,1);
    for i = 1:cd_tot
        basins(i) = test_basin(lawbook(i,:));
    end
    M_temp = max(M_temp,max(basins));

    code_swapping = zeros(cd_tot,max(faithful_set));
    for i = 1:cd_tot
        code_swapping(i,lawbook(i,:)) = true;
    end
    code_swapping = 1-code_swapping;

    hold on

    imagesc(ax(k),1:cd_tot,1:max(faithful_set),code_swapping.');
    colormap([0 0 0; 1 1 1]);

    x = (1:cd_tot);
    yyaxis(ax(k),'right')
    plot(ax(k),x,basins,'r-o','LineWidth',2)

    title(sprintf("Codebooks by %s",algo_titles(k)))

end

for k = 1:4
    yyaxis(ax(k),'right')
    ax(k).YColor = 'k';
    ylim(ax(k),[0,M_temp+1])
    xlim(ax(k),[0,cd_tot+1])
end

t.OuterPosition = [0.05 0.05 0.85 0.9];

linkaxes(ax,'x')
title(t,"Generation of Codebooks")
xlabel(t,"Codebook #")
ylabel(t,"Code (Walsh) #")

annotation('textbox',[0.88 0.57 0.1 0.1],...
    'String','Min Basin Size',...
    'Rotation',-90,...
    'EdgeColor','none',...
    'HorizontalAlignment','center', ...
    'FontSize',12);

saveas(t,'comparing_algos.png')

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
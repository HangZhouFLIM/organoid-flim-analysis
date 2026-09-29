function crypt_to_villus_gradients_manual_1
% Interactive pipeline for crypt->villus intensity gradients
% - Select multiple CSV files (each with Id, MeanInt, CentroidX, CentroidY)
% - For each file, enter crypt cell ID and villus cell ID
% - Computes projection, bins, summary stats
% - Produces & saves plots (Prism-like styling):
%     * Mean + SEM (pooled)
%     * Combined single-cell heatmap (all samples)
%     * Per-sample panel heatmap (one row per sample)
%     * Bin-averaged intensity heatmap
%     * Boxplots per bin
% - Exports CSV files for each plot
%
% Output folder: ./output
%

%% --------------------- USER SETTINGS ---------------------
numBins        = 50;
pctList        = [10 25 50 75 90];
outFolder      = './output';
perSampleAlpha = 0.25; %#ok<NASGU> % (reserved; not used in this version)
lowessSpan     = 0.2;

if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% --------------------- SELECT FILES + AXES ---------------------
[fileList, axisDef] = chooseFilesAndAxes();
if isempty(fileList), disp('No files selected. Exiting.'); return; end

%% --------------------- PROCESS SAMPLES ---------------------
allRows = table();
perSampleSummary = table();
sampleInfo = struct();

for k = 1:numel(fileList)
    fpath = fileList{k};
    [~,fname,ext] = fileparts(fpath);
    sampleName = [fname ext];
    T = readtable(fpath);

    required = {'Id','MeanInt','CentroidX','CentroidY'};
    assert(all(ismember(required, T.Properties.VariableNames)), ...
        'File %s must contain columns: %s', fpath, strjoin(required, ', '));

    cryptId  = axisDef{k}.cryptId;
    villusId = axisDef{k}.villusId;
    [baseXY, tipXY] = getAxisFromIds(T, cryptId, villusId);

    v = (tipXY - baseXY);
    vlen2 = sum(v.^2);
    assert(vlen2 > 0, 'Chosen crypt and villus points are identical for %s', sampleName);

    pts = [T.CentroidX, T.CentroidY];
    proj = ((pts - baseXY) * v.') / vlen2;
    proj01 = min(max(proj,0),1);

    edges = linspace(0,1,numBins+1);
    [~,~,binIdx] = histcounts(proj01, edges);
    binCenters = edges(1:end-1) + diff(edges)/2;

    T.Sample = repmat(string(sampleName), height(T), 1);
    T.Proj01 = proj01;
    T.Bin    = binIdx;

    allRows = [allRows; T]; %#ok<AGROW>

    S = summarizeBins(T, binCenters, pctList);
    S.Sample = repmat(string(sampleName), height(S), 1);
    perSampleSummary = [perSampleSummary; S]; %#ok<AGROW>

    [~,ord] = sort(T.Proj01);
    sampleInfo(k).sampleName = sampleName; 
    sampleInfo(k).Proj01 = T.Proj01(ord);
    sampleInfo(k).MeanInt = T.MeanInt(ord);

    sampleTableOrdered = table(T.Id(ord), T.Proj01(ord), T.MeanInt(ord), ...
        'VariableNames', {'Id','Proj01','MeanInt'});
    writetable(sampleTableOrdered, fullfile(outFolder, sprintf('heatmap_%s_ordered.csv', matlab.lang.makeValidName(sampleName))));
end

writetable(allRows, fullfile(outFolder,'all_cells_with_projection.csv'));

%% --------------------- OVERALL SUMMARY ---------------------
edges = linspace(0,1,numBins+1);
binCenters = edges(1:end-1) + diff(edges)/2;
overall = summarizeBins(allRows, binCenters, pctList);
writetable(overall, fullfile(outFolder,'bin_summary_overall.csv'));
writetable(perSampleSummary, fullfile(outFolder,'bin_summary_per_sample.csv'));

boxLong = allRows(:, {'Bin','Sample','Id','MeanInt'});
writetable(boxLong, fullfile(outFolder,'boxplot_long_per_bin.csv'));

[~, ordAll] = sort(allRows.Proj01);
heat_all = table(allRows.Sample(ordAll), allRows.Id(ordAll), allRows.Proj01(ordAll), allRows.MeanInt(ordAll), ...
    'VariableNames', {'Sample','Id','Proj01','MeanInt'});
writetable(heat_all, fullfile(outFolder,'heatmap_allcells_ordered.csv'));

lowessX = overall.BinCenter;
lowessY = overall.Median;
smY = smooth(lowessX, lowessY, lowessSpan, 'lowess');
lowessT = table(lowessX, smY, 'VariableNames', {'BinCenter','SmoothedMedian'});
writetable(lowessT, fullfile(outFolder,'lowess_smoothed.csv'));

%% --------------------- PLOTS (Prism-like) ---------------------
makePlotsAndHeatmaps(allRows, overall, sampleInfo, numBins, outFolder);

end

%% --------------------- INTERACTIVE FILE + AXIS SELECTION ---------------------
function [files, axisDef] = chooseFilesAndAxes()
    [fn,pth] = uigetfile('*.csv','Select CSV files','MultiSelect','on');
    if isequal(fn,0), files={}; axisDef={}; return; end
    if ischar(fn), fn={fn}; end
    files = fullfile(pth, fn);
    axisDef = cell(numel(files),1);
    for k=1:numel(files)
        T = readtable(files{k});
        [~,basename,ext] = fileparts(files{k});
        sampName = [basename ext];
        fprintf('\n--- File %d: %s ---\n', k, sampName);
        previewTbl = T(1:min(height(T),30), {'Id','CentroidX','CentroidY'});
        disp(previewTbl);
        cryptId = input(sprintf('Enter crypt_base Id for %s: ', sampName));
        villusId = input(sprintf('Enter villus_tip Id for %s: ', sampName));
        if ~any(T.Id==cryptId), error('cryptId %d not found', cryptId); end
        if ~any(T.Id==villusId), error('villusId %d not found', villusId); end
        centroids=[T.CentroidX,T.CentroidY];
        f = figure('Name',sprintf('Centroids preview: %s',sampName),'Color','w');
        styleDefaults(f); % apply default fonts etc. to preview too
        scatter(centroids(:,1),centroids(:,2),12,'k','filled'); hold on;
        rC=T(T.Id==cryptId,:); rV=T(T.Id==villusId,:);
        scatter(rC.CentroidX,rC.CentroidY,80,'b','filled','MarkerEdgeColor','w');
        scatter(rV.CentroidX,rV.CentroidY,80,'r','filled','MarkerEdgeColor','w');
        legend({'cells','crypt','villus'},'Box','off','Location','best');
        axis equal; set(gca,'YDir','reverse');
        title(sampName,'FontWeight','bold');
        prismAxes(gca);
        input('Press Enter to continue...','s'); close(f);
        axisDef{k}=struct('cryptId',cryptId,'villusId',villusId);
    end
end

%% --------------------- AXIS LOOKUP ---------------------
function [baseXY, tipXY] = getAxisFromIds(T, cryptId, villusId)
    rowC = T(T.Id==cryptId,:); rowV = T(T.Id==villusId,:);
    baseXY=[rowC.CentroidX,rowC.CentroidY];
    tipXY=[rowV.CentroidX,rowV.CentroidY];
end

%% --------------------- BIN SUMMARY ---------------------
function S = summarizeBins(T, binCenters, pctList)
    numBins = numel(binCenters);
    Mean=nan(numBins,1); SEM=nan(numBins,1); N=zeros(numBins,1);
    Pcts=nan(numBins,numel(pctList));
    for b=1:numBins
        mask=(T.Bin==b); x=T.MeanInt(mask); x=x(~isnan(x));
        N(b)=numel(x);
        if N(b)>0
            Mean(b)=mean(x); SEM(b)=std(x)/sqrt(N(b));
            Pcts(b,:)=prctile(x,pctList);
        end
    end
    S=table((1:numBins).',binCenters(:),N(:),Mean(:),SEM(:), ...
        'VariableNames',{'Bin','BinCenter','N','Mean','SEM'});
    for i=1:numel(pctList)
        S.(sprintf('P%02d',pctList(i)))=Pcts(:,i);
    end
    if any(strcmp('P50',S.Properties.VariableNames))
        S.Median=S.P50; else S.Median=NaN(size(Mean)); end
end

%% --------------------- PLOTTING ---------------------
function makePlotsAndHeatmaps(allRows, overall, sampleInfo, numBins, outFolder)
    % Common color limits for intensity-based plots
    [~,ordAll]=sort(allRows.Proj01); valsAll=allRows.MeanInt(ordAll);
    clim = prctile(valsAll,[5 95]);

    %% 1) Mean + SEM (pooled; Prism-like)
    f1 = figure('Color','w','Position',[200 200 1000 600]); styleDefaults(f1);
    ax1 = axes('Parent',f1); hold(ax1,'on'); box(ax1,'off');
    p = plot(ax1, overall.BinCenter, overall.Mean, 'k-', 'LineWidth', 2);
    fb = fillBetween(ax1, overall.BinCenter, overall.Mean-overall.SEM, overall.Mean+overall.SEM, [0 0 0], 0.15);
    uistack(fb,'bottom'); uistack(p,'top');
    xlabel(ax1,'Normalized position (0 = crypt, 1 = villus)','FontWeight','bold');
    ylabel(ax1,'Mean intensity','FontWeight','bold');
    title(ax1,'Pooled mean ± SEM','FontWeight','bold');
    prismAxes(ax1);
    export600(f1, fullfile(outFolder,'fig_mean_sem.png'));

    %% 2) Combined single-cell heatmap
    f2 = figure('Color','w','Position',[200 200 1200 320]); styleDefaults(f2);
    ax2 = axes('Parent',f2); hold(ax2,'on'); box(ax2,'off');
    H = repmat(valsAll.',24,1);
    imagesc(ax2, [0 1],[0 1], H);
    set(ax2,'YTick',[],'YColor','none'); % de-clutter
    xlabel(ax2,'Normalized position (all samples)','FontWeight','bold');
    title(ax2,'Single-cell intensities (combined)','FontWeight','bold');
    c = colorbar(ax2); ylabel(c,'Mean intensity','FontWeight','bold');
    caxis(ax2,clim); colormap(ax2,parula);
    prismAxes(ax2,true); % heatmap-friendly adjustments
    export600(f2, fullfile(outFolder,'fig_heatmap_allcells.png'));

    %% 3) Per-sample panel heatmap
    ns = numel(sampleInfo);
    f3 = figure('Color','w','Position',[200 200 1200, max(300, 160+130*ns)]); styleDefaults(f3);
    for i=1:ns
        ax = subplot(ns,1,i,'Parent',f3);
        prismAxes(ax,true);
        vals_i = sampleInfo(i).MeanInt;
        if isempty(vals_i)
            imagesc(ax, [0 1],[0 1], zeros(10,1));
        else
            imagesc(ax, [0 1],[0 1], repmat(vals_i.',20,1));
        end
        set(ax,'YTick',[],'YColor','none');
        title(ax, sampleInfo(i).sampleName, 'Interpreter','none','FontWeight','bold');
        colormap(ax,parula); caxis(ax,clim);
        if i==ns
            xlabel(ax,'Normalized position (within sample)','FontWeight','bold');
        else
            set(ax,'XTickLabel',[]);
        end
    end
    % Safe group title with fallback if sgtitle not available
    try
        sgtitle(f3,'Per-sample single-cell heatmaps','FontWeight','bold');
    catch
        % no-op if sgtitle unavailable
    end
    export600(f3, fullfile(outFolder,'fig_heatmap_per_sample_panel.png'));

    %% 4) Bin-averaged intensity heatmap
    f4 = figure('Color','w','Position',[200 200 1200 320]); styleDefaults(f4);
    ax4 = axes('Parent',f4); prismAxes(ax4,true); hold(ax4,'on');
    binMeans = overall.Mean(:)';
    H2 = repmat(binMeans,20,1);
    imagesc(ax4, [0 1],[0 1], H2);
    xlabel(ax4,'Normalized position (bin centers)','FontWeight','bold');
    set(ax4,'YTick',[],'YColor','none');
    title(ax4,'Bin-averaged intensities heatmap','FontWeight','bold');
    c = colorbar(ax4); ylabel(c,'Mean intensity','FontWeight','bold');
    caxis(ax4,clim); colormap(ax4,parula);
    export600(f4, fullfile(outFolder,'fig_heatmap_bin_means.png'));

    %% 5) Boxplot per bin
    f5 = figure('Color','w','Position',[200 200 1400 520]); styleDefaults(f5);
    ax5 = axes('Parent',f5); hold(ax5,'on'); box(ax5,'off');
    valid = ~isnan(allRows.Bin) & allRows.Bin > 0;
    boxplot(ax5, allRows.MeanInt(valid), allRows.Bin(valid), ...
        'PlotStyle','compact', 'Labels',[]);
    % Style boxplot lines a bit bolder
    set(findobj(ax5,'Tag','Box'),'LineWidth',1.2);
    set(findobj(ax5,'Tag','Median'),'LineWidth',1.2);
    set(findobj(ax5,'Tag','Whisker'),'LineWidth',1.2);
    set(findobj(ax5,'Tag','Outliers'),'MarkerSize',3);

    edges = linspace(0,1,numBins+1);
    xticksLoc = round(linspace(1,numBins,11));
    xticks(ax5, xticksLoc);
    xticklabels(ax5, arrayfun(@(x)sprintf('%.2f',edges(x)), xticksLoc,'UniformOutput',false));
    xlabel(ax5,'Normalized position (bin index/edges)','FontWeight','bold');
    ylabel(ax5,'Intensity','FontWeight','bold');
    title(ax5,'Per-bin intensity distributions','FontWeight','bold');
    prismAxes(ax5);
    export600(f5, fullfile(outFolder,'fig_boxplots_per_bin.png'));
end

%% --------------------- STYLING HELPERS ---------------------
function styleDefaults(figHandle)
    % Per-figure defaults (avoids altering global graphics root)
    set(figHandle, ...
        'DefaultAxesFontName','Arial', ...
        'DefaultAxesFontSize',12, ...
        'DefaultAxesFontWeight','bold', ...
        'DefaultAxesLineWidth',2.5, ...
        'DefaultAxesTickDir','out', ...
        'DefaultAxesTickLength',[0.015 0.015], ...
        'DefaultAxesBox','off', ...
        'DefaultLineLineWidth',2, ...
        'DefaultFigurePaperPositionMode','auto');
end

function prismAxes(ax, isHeatmap)
    if nargin<2, isHeatmap=false; end
    set(ax, 'Box','off', ...              % removes top/right frame
            'TickDir','out', ...
            'LineWidth',2.5, ...
            'FontWeight','bold', ...
            'FontName','Arial', ...
            'FontSize',12, ...
            'Layer','top');               % ticks above images
    if isHeatmap
        % Keep axes clean for heatmaps
        grid(ax,'off');
        % Thin axis lines look cleaner under images
        ax.LineWidth = 1.0;
    end
end

function export600(figHandle, outPath)
    % Save at 600 dpi. Try exportgraphics (newer), fall back to print (older).
    [~,~,ext] = fileparts(outPath);
    if isempty(ext), outPath = [outPath '.png']; end
    try
        exportgraphics(figHandle, outPath, 'Resolution', 600);
    catch
        print(figHandle, outPath, '-dpng', '-r600');
    end
end

function h=fillBetween(ax,x,y1,y2,colorRGB,alphaVal)
    % Fill between y1 and y2 over x on axes ax
    x=x(:)'; y1=y1(:)'; y2=y2(:)';
    hold(ax,'on');
    h=fill(ax, [x fliplr(x)], [y1 fliplr(y2)], colorRGB, ...
        'EdgeColor','none', 'FaceAlpha', alphaVal);
end

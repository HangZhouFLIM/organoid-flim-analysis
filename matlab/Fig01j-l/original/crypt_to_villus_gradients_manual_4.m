function crypt_to_villus_gradients_manual_1
% Interactive pipeline for crypt->villus intensity gradients
% - Uses column 'MeanTauModulation'
% - Produces & saves plots (Prism-like):
%     * Mean + SEM (pooled)
%     * Combined single-cell heatmap
%
% Output folder: ./output

%% --------------------- USER SETTINGS ---------------------
numBins        = 50;
pctList        = [10 25 50 75 90];
outFolder      = './output';
lowessSpan     = 0.2;

if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% --------------------- SELECT FILES + AXES ---------------------
[fileList, axisDef] = chooseFilesAndAxes();
if isempty(fileList), disp('No files selected. Exiting.'); return; end

%% --------------------- PROCESS SAMPLES ---------------------
allRows = table();
sampleInfo = struct();

for k = 1:numel(fileList)
    fpath = fileList{k};
    [~,fname,ext] = fileparts(fpath);
    sampleName = [fname ext];
    T = readtable(fpath);

    required = {'Id','MeanTauModulation','CentroidX','CentroidY'};
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

    [~,ord] = sort(T.Proj01);
    sampleInfo(k).sampleName = sampleName; %#ok<SAGROW>
    sampleInfo(k).Proj01 = T.Proj01(ord);
    sampleInfo(k).MeanTauModulation = T.MeanTauModulation(ord);

    sampleTableOrdered = table(T.Id(ord), T.Proj01(ord), T.MeanTauModulation(ord), ...
        'VariableNames', {'Id','Proj01','MeanTauModulation'});
    writetable(sampleTableOrdered, fullfile(outFolder, sprintf('heatmap_%s_ordered.csv', matlab.lang.makeValidName(sampleName))));
end

writetable(allRows, fullfile(outFolder,'all_cells_with_projection.csv'));

%% --------------------- OVERALL SUMMARY ---------------------
binCenters = linspace(0,1,numBins+1);
binCenters = binCenters(1:end-1) + diff(binCenters)/2;
overall = summarizeBins(allRows, binCenters, pctList);
writetable(overall, fullfile(outFolder,'bin_summary_overall.csv'));

[~, ordAll] = sort(allRows.Proj01);
heat_all = table(allRows.Sample(ordAll), allRows.Id(ordAll), allRows.Proj01(ordAll), allRows.MeanTauModulation(ordAll), ...
    'VariableNames', {'Sample','Id','Proj01','MeanTauModulation'});
writetable(heat_all, fullfile(outFolder,'heatmap_allcells_ordered.csv'));

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
        styleDefaults(f);
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
        mask=(T.Bin==b); x=T.MeanTauModulation(mask); x=x(~isnan(x));
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
    [~,ordAll]=sort(allRows.Proj01); valsAll=allRows.MeanTauModulation(ordAll);
    clim = prctile(valsAll,[5 95]);

    %% 1) Mean + SEM (pooled; Prism-like)
    f1 = figure('Color','w','Position',[200 200 1000 600]); styleDefaults(f1);
    ax1 = axes('Parent',f1); hold(ax1,'on'); box(ax1,'off');
    plot(ax1, overall.BinCenter, overall.Mean, 'k-', 'LineWidth', 2);
    fillBetween(ax1, overall.BinCenter, overall.Mean-overall.SEM, overall.Mean+overall.SEM, [0 0 0], 0.15);
    xlabel(ax1,'Normalized position (0 = crypt, 1 = villus)','FontWeight','bold');
    ylabel(ax1,'MeanTauModulation','FontWeight','bold');
    title(ax1,'Pooled mean ± SEM','FontWeight','bold');
    prismAxes(ax1);
    export600(f1, fullfile(outFolder,'fig_mean_sem.png'));

    %% 2) Combined single-cell heatmap
    f2 = figure('Color','w','Position',[200 200 1200 320]); styleDefaults(f2);
    ax2 = axes('Parent',f2); hold(ax2,'on'); box(ax2,'off');
    H = repmat(valsAll.',24,1);
    imagesc(ax2, [0 1],[0 1], H);
    set(ax2,'YTick',[],'YColor','none');
    xlabel(ax2,'Normalized position (all samples)','FontWeight','bold');
    title(ax2,'Single-cell MeanTauModulation (combined)','FontWeight','bold');
    c = colorbar(ax2); ylabel(c,'MeanTauModulation','FontWeight','bold');
    caxis(ax2,clim); colormap(ax2,parula);
    prismAxes(ax2,true);
    export600(f2, fullfile(outFolder,'fig_heatmap_allcells.png'));
end

%% --------------------- STYLING HELPERS ---------------------
function styleDefaults(figHandle)
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
    set(ax, 'Box','off','TickDir','out','LineWidth',2.5,'FontWeight','bold','FontName','Arial','FontSize',12,'Layer','top');
    if isHeatmap
        grid(ax,'off'); ax.LineWidth = 1.0;
    end
end

function export600(figHandle, outPath)
    [~,~,ext] = fileparts(outPath);
    if isempty(ext), outPath = [outPath '.png']; end
    try
        exportgraphics(figHandle, outPath, 'Resolution', 600);
    catch
        print(figHandle, outPath, '-dpng', '-r600');
    end
end

function h=fillBetween(ax,x,y1,y2,colorRGB,alphaVal)
    x=x(:)'; y1=y1(:)'; y2=y2(:)';
    hold(ax,'on');
    h=fill(ax, [x fliplr(x)], [y1 fliplr(y2)], colorRGB, 'EdgeColor','none', 'FaceAlpha', alphaVal);
end

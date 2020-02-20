%% Generate acq name
BaseStr = regexprep([char(ispc.*'Z:\Images2019\') char(isunix.*'/bigstore/Images2019/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'TNFTitr_HighMOI_Dec122019_2019Dec12';
acquisition = 2;
%% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath,[],1);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));

%%
   set(0,'DefaultTextInterpreter', 'tex')
set(0, 'DefaultAxesFontName', 'Arial')
set(0, 'DefaultUIControlFontName', 'Arial')
set(0,'defaulttextfontname','Arial');
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesColor','w')
set(groot,'defaultAxesTickLength',[0.03 0.01])
set(groot,'defaultLineLineWidth',2)
%% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
dt = 1/3;%h

cTNF = ([100*(1/2).^[0:8], 0]);
%% Load Results

R = MultiPositionSingleCellVirusResults(fpath)
Wells = R.PosNames;
frames = R.Frames;
t = [R.Frames-1]'.*dt;

%% Calibrate TNF to infected death rates
allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})
allParams = allParams(:,2);
allParams = reshape(allParams,[],4);

close all
figure('color','w','Position',[100,100, 300, 300])
axes('Position', [0.18, 0.14, 0.73, 0.75])
x = cTNF'
y = mean(allParams(:,1:2),2);
sigmay = std(allParams(:,1:2),[],2);


h = ploterr(cTNF,y,[], sigmay,'ko', 'logx');
for jj=1:numel(h)
h(jj).LineWidth = 2;
end
%sigmay = std(allParams(:,1:2),[],2)./y;
%sigmay = ones(size(y));
BETA0 = [max(y), mean(x), 0];
BETA = lsqcurvefit(@(beta, xdata) HillFunction(beta, xdata), BETA0, x, y);
hold on
x = 0:0.1:max(x)+50;
y = HillFunction(BETA,x);
h = plot(x,y);
h.LineWidth = 2;

set(gca,'XLim',[0.2 200])

shg
HillFuncBeta = BETA;
ylabel('Infected death rate - \beta_{{i}}')
xlabel('[TNF] ng/ml')
%R.setData('HillFuncBeta', HillFuncBeta)


%%
savpath = '/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/ModelFigures/'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-deps','-r600',[savpath 'DeathVsTNFHill']);

%% Plot S/A tradeoff
close all
figure('color','w','Position',[100,100, 300, 300])
axes('Position', [0.18, 0.14, 0.73, 0.75])
allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})

allParams = allParams(:,2);
allParams = reshape(allParams,[],4);

tzeva = flipud(parula(numel(cTNF)));

for i=1:numel(cTNF)
h = ploterr(mean(allParams(i,1:2),2), mean(allParams(i,3:4),2), std(allParams(i,1:2),[],2), std(allParams(i,3:4),[],2),'o');
hold on
for jj=1:numel(h)
h(jj).Color = tzeva(i,:);
h(jj).MarkerFaceColor = tzeva(i,:);
h(jj).LineWidth = 2;
end

end
set(gca,'YLim',[-0.001, 0.02],'XLim',[-0.005, 0.12], 'YTick', [0 0.0050 0.0100 0.0150 0.0200], 'YTickLabel', num2str([0 0.0050 0.0100 0.0150 0.0200]','%.3G'))
xlabel('Infected death rate - \beta_{{i}}')
ylabel('Bystander death rate - \beta_{{b}}')


%% Calibrate S/A tradeoff
tzeva = flipud(parula(numel(cTNF)));

for i=1:numel(cTNF)
h = ploterr(mean(allParams(i,1:2),2), mean(allParams(i,3:4),2), std(allParams(i,1:2),[],2), std(allParams(i,3:4),[],2),'o');
hold on
for jj=1:numel(h)
h(jj).Color = tzeva(i,:);
h(jj).MarkerFaceColor = tzeva(i,:);
h(jj).LineWidth = 2;
end

end
BETA0 = [0.2, 0.01];
%BETA = lsqcurvefit(@Rectifier, BETA0, mean(allParams(:,1:2),2), mean(allParams(:,3:4),2));
hold on
%x=0:0.001:0.13;
%plot(x, Rectifier(BETA,x))
%shg
x = mean(allParams(:,1:2),2);
y = mean(allParams(:,3:4),2);
sigmax = std(allParams(:,1:2),[],2);
sigmay = std(allParams(:,3:4),[],2);

BETA = lsqcurvefit(@(beta, xdata) Rectifier(beta, xdata)./sigmay, BETA0, x, y./sigmay)
x=0:0.001:0.1;
y=Rectifier(BETA,x);
h = plot(x, y,'-r');
xlabel('Infected death rate - \beta_{{i}}')
ylabel('Bystander death rate - \beta_{{b}}')
shg

RectifierBETA = BETA;
R.setData('RectifierBETA', RectifierBETA)

%%
savpath = '/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/ModelFigures/'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath 'SAToff']);

%% This is how you calculate rates:
HillFuncBeta = R.getData('HillFuncBeta')
RectifierBETA = R.getData('RectifierBETA')
TNF = 100; %input TNF
infDeathRate = HillFunction(HillFuncBeta,TNF);
basalDeathRate = Rectifier(RectifierBETA,infDeathRate);
%% Make triangular grid
A = triangleGrid([-10 -10 10 10], [0,0], 1);
nCells = size(A,1);
%find all nearest neighnors and create sparse matrix
[IDX, D] = rangesearch(A, A,1.01);
IDX = cellfun(@(x,y) x(y>0), IDX,D,'uniformoutput',false);
pointsToMatch = arrayfun(@(x,y) repmat(x,y,1),1:numel(IDX), cellfun(@numel,IDX)','UniformOutput', false);
nnMatrix = sparse(cat(1,pointsToMatch{:}), cat(2,IDX{:})',1);

%% Init cells, small frac infected
%fracInf = 0.01;
%infected = rand(nCells,1)<fracInf;
infected = zeros(nCells,1);
infected(116)=1;

x= [~infected, infected, zeros(nCells,1), zeros(nCells,1)]';
x = x(:);

figure('color','w')
ax1 = axes('position',[0,0,1,1]);
%voronoi(A(:,1), A(:,2)); hold on;
h(1) = scatter(A(logical(x(1:4:end)),1),A(logical(x(1:4:end)),2),100,'b','filled')
hold on;
h(2) = scatter(A(logical(x(2:4:end)),1),A(logical(x(2:4:end)),2),100,'r','filled')
h(3) = scatter(A(logical(x(3:4:end)),1),A(logical(x(3:4:end)),2),100,'k','filled')
h(4) = scatter(A(logical(x(4:4:end)),1),A(logical(x(4:4:end)),2),100,'g','filled')
axis equal
set(ax1,'xcolor','w','ycolor','w','xlim', [-10, 10],'ylim', [-10, 10])
shg
%hl = legend(h,'Healthy','Infected','False positive','Dead post infection')
%hl.Location = 'northeastoutside';
%hl.Box = 'off';

figure(gcf)
%set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
%print(gcf,'-dpng','-r300',['/Users/Alonyan/Desktop/' 'Initial']);

%%
savpath = '/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/ModelFigures/'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath 'SAToff']);

%% Run model

VGR = 1/10;%viral growth rate
VI = 1/4; %viral infectivity

TNF = 0;
infDeathRate = HillFunction(HillFuncBeta,TNF)
basalDeathRate = Rectifier(RectifierBETA,infDeathRate)


xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2))

%
figure('color','w')
ax1 = axes('position',[0,0,1,1]);
scatter(A(logical(xOut(1:4:end)),1),A(logical(xOut(1:4:end)),2),100,'b','filled');
hold on;
scatter(A(logical(xOut(2:4:end)),1),A(logical(xOut(2:4:end)),2),100,'r','filled');
scatter(A(logical(xOut(3:4:end)),1),A(logical(xOut(3:4:end)),2),100,'k','filled');
scatter(A(logical(xOut(4:4:end)),1),A(logical(xOut(4:4:end)),2),100,'g','filled');
axis equal
set(ax1,'xcolor','w','ycolor','w');
an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['TNF=' num2str(TNF) 'ng/ml']},'LineStyle','none');

shg
hl = legend('Healthy','Infected','False positive','Dead post infection');
hl.Location = 'northeastoutside';
hl.Box = 'off';
%%
figure(gcf)
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-dpng','-r300',['/Users/Alonyan/Desktop/' 'HighAlpha2']);








%% Make triangular grid
A = triangleGrid([-10 -10 10 10], [0,0], 1);
nCells = size(A,1);
%find all nearest neighnors and create sparse matrix
[IDX, D] = rangesearch(A, A,1.01);
IDX = cellfun(@(x,y) x(y>0), IDX,D,'uniformoutput',false);
pointsToMatch = arrayfun(@(x,y) repmat(x,y,1),1:numel(IDX), cellfun(@numel,IDX)','UniformOutput', false);
nnMatrix = sparse(cat(1,pointsToMatch{:}), cat(2,IDX{:})',1);

%% Get Statistics
VGR = 1/10;%viral growth rate
VI = 1/4; %viral infectivity

nRuns = 500;
TNF = [0 fliplr(100*(1/sqrt(2)).^[0:17])];
Stats = {}
for j=1:numel(TNF)
    infDeathRate = HillFunction(HillFuncBeta,TNF(j))
    basalDeathRate = Rectifier(RectifierBETA,infDeathRate)
    pStats = cell(nRuns,1);
    parfor i=1:nRuns
        pStats{i} = modelRunForStats(nnMatrix,[],VGR,VI,basalDeathRate,infDeathRate);
    end
    Stats{j} = cat(2,pStats{:})';
end
%%
R.setData('modelTNF', TNF)
R.setData('modelSTATs', Stats)
R.saveResults


%% Plot Model statistics
close all
a = cellfun(@(x) 100*x(:,1), Stats,'uniformOutput',false)
yyaxis left
plot(median(cat(2,a{:})),'b');
set(gca,'YLim',[0,100])
ylabel('%Healthy at 72h')

b = cellfun(@(x) 100*x(:,2), Stats,'uniformOutput',false)
yyaxis right
plot(median(cat(2,b{:})),'r');
set(gca,'YLim',[0,100])
ylabel('%Infected at 72h')
hold on
c = cellfun(@(x) 100*x(:,4), Stats,'uniformOutput',false)
plot(median(cat(2,c{:})),'g');
set(gca,'YLim',[0,100])
c = cellfun(@(x) 100*x(:,3), Stats,'uniformOutput',false)
plot(median(cat(2,c{:})),'k');
set(gca,'YLim',[0,100])

hold on

h = boxplot(cat(2,a{:}), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', 'b');
shg
xlabel('[TNF] ng/ml')
set(gca,'YLim',[-1,100])
h = boxplot(cat(2,b{:}), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', 'r');
shg
set(gca,'YLim',[-1,100])

%%
binEdges = 0.05:0.05:1;
figure
subplot(2,2,1)
histogram(StatsNoTNF(:,1),'BinEdges',binEdges)
hold on
histogram(StatsHighTNF(:,1),'BinEdges',binEdges)
%histogram(StatsTooMuchTNF(:,1),'BinEdges',binEdges)
title('Healthy')

subplot(2,2,2)
histogram(StatsNoTNF(:,2),'BinEdges',binEdges)
hold on
histogram(StatsHighTNF(:,2),'BinEdges',binEdges)
%histogram(StatsTooMuchTNF(:,2),'BinEdges',binEdges)
title('Infected @ 72h')

subplot(2,2,3)
histogram(StatsNoTNF(:,3),'BinEdges',binEdges)
hold on
histogram(StatsHighTNF(:,3),'BinEdges',binEdges)
%histogram(StatsTooMuchTNF(:,3),'BinEdges',binEdges)
title('Uninfected Dead')

subplot(2,2,4)
histogram(StatsNoTNF(:,4),'BinEdges',binEdges)
hold on
histogram(StatsHighTNF(:,4),'BinEdges',binEdges)
%histogram(StatsTooMuchTNF(:,4),'BinEdges',binEdges)
title('Infected Dead')







%% MakeFig for grant

   set(0,'DefaultTextInterpreter', 'tex')
set(0, 'DefaultAxesFontName', 'Arial')
set(0, 'DefaultUIControlFontName', 'Arial')
set(0,'defaulttextfontname','Arial');

fracInf = 0.005;
infected = rand(nCells,1)<fracInf;

x= [~infected, infected, zeros(nCells,1), zeros(nCells,1)]';
x = x(:);

figure('Units' ,'Inches','Position',[1, 1 4 2],'color','w')


VGR = 1/6;%viral growth rate
VI = 1/2; %viral infectivity

basalDeathRate=1/1000;%effect of TNF on deaths
infDeathRate = 1/24;

xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2))

ax2 = axes('position',[0.22,0.01,0.22,1])

h(1) = scatter(A(logical(xOut(1:4:end)),1),A(logical(xOut(1:4:end)),2),6,'b','filled');
hold on;
h(2) = scatter(A(logical(xOut(2:4:end)),1),A(logical(xOut(2:4:end)),2),6,'r','filled');
h(3) = scatter(A(logical(xOut(3:4:end)),1),A(logical(xOut(3:4:end)),2),6,[0.5 0.5 0.5],'filled');
h(4) = scatter(A(logical(xOut(4:4:end)),1),A(logical(xOut(4:4:end)),2),6,'g','filled');
axis equal
set(ax2,'xcolor','none','ycolor','none');
tl = title('No TNF')
tl.Position=[0.48, 11, 0];

%an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['VI=' num2str(VI)], ['\alpha=' num2str(alpha)], ['VGR=' num2str(VGR)]},'LineStyle','none');

shg
% hl.Location = 'northeastoutside';
% hl.Box = 'off';
hl = legend(h(1:2),' Healthy',' Infected');
hl.Location = 'northeastoutside';
hl.Box = 'off';
hl.FontSize = 8;
hl.Position = [0.24, 0.1, 0.1, 0.1]



ax3 = axes('position',[0.46,0.01,0.22,1])


basalDeathRate=1/100;%effect of TNF on deaths
infDeathRate = 1/4;

xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2))


h(1) = scatter(A(logical(xOut(1:4:end)),1),A(logical(xOut(1:4:end)),2),6,'b','filled');
hold on;
h(2) = scatter(A(logical(xOut(2:4:end)),1),A(logical(xOut(2:4:end)),2),6,'r','filled');
h(3) = scatter(A(logical(xOut(3:4:end)),1),A(logical(xOut(3:4:end)),2),6,[0.5 0.5 0.5],'filled');
h(4) = scatter(A(logical(xOut(4:4:end)),1),A(logical(xOut(4:4:end)),2),6,'g','filled');
axis equal
set(ax3,'xcolor','none','ycolor','none');
%an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['VI=' num2str(VI)], ['\alpha=' num2str(alpha)], ['VGR=' num2str(VGR)]},'LineStyle','none');
tl = title('TNF');
tl.Position=[0.5, 11, 0];
shg
hl = legend(h(3:4),' Uninfected-Dead',' Infected-Dead');
hl.Location = 'northeastoutside';
hl.Box = 'off';
hl.FontSize = 8;
hl.Position = [0.5, 0.1, 0.1, 0.1]



ax4 = axes('position',[0.79, 0.28 0.178, 0.45])
binEdges = -0.025:0.05:1.025;

histogram(StatsNoTNF(:,1),'BinEdges',binEdges)
hold on
histogram(StatsHighTNF(:,1),'BinEdges',binEdges)
set(gca,'ytick', 0:100:500, 'ytickLabel',100*[0:100:500]./nRuns,'xtick',[0:0.25:1],'xticklabel',100*[0:0.25:1],'xlim',[-0.03, 1],'fontsize',8,'ticklength',[0.02, 0.1])
xlabel('% healthy at 72h','fontsize',8)
ylabel('% of simulations','fontsize',8)
[hl, ~]= legend({'No TNF', 'TNF'})
hl.Box = 'off';
hl.FontSize = 8;
hl.Position = [0.82, 0.81, 0.1 0.1]
tzeva = lines(2);
annotation('rectangle',[0.83,0.88,0.05,0.05],'facecolor',tzeva(1,:)','facealpha',0.6)
annotation('rectangle',[0.83,0.8,0.05,0.05],'facecolor',tzeva(2,:)','facealpha',0.6)

%hLine = findobj(objh,'type','patch')
%        set(hLine,'XData', [0.27 0.27 0.35 0.35 0.27]);


figure(gcf)


ImSketch = imread('/bigstore/GeneralStorage/Alon/Figures/R01June2019/Fig2Cartoon.tiff');

width = 0.4;
   axes('Position',[-0.1 0.27 width width.*(size(ImSketch,1)./size(ImSketch,2))])

h= imshow(ImSketch(:,:,1:3)); figure(gcf)


   a1=annotation('textbox',[0 0.95 0.05 0.05],'string','a', 'Fontsize',10,'LineStyle','none')
   a2 = annotation('textbox',[0.2 0.95 0.05 0.05],'string','b', 'Fontsize',10,'LineStyle','none')
   
   a3=annotation('textbox',[ 0.68 0.95    0.05    0.05],'string','c', 'Fontsize',10,'LineStyle','none')
   
   
%%
figure(gcf)
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','Renderer','opengl')
print(gcf,'-depsc2','-r600',['/bigstore/GeneralStorage/Alon/Figures/R01June2019/' 'Fig2']);
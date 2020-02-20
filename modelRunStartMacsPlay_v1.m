
load('/bigstore/GeneralStorage/Alon/Figures/R01June2019/SingleTP.mat')
%% Make spatial arangement of cells
A = SingleTP.Centroids;%(1000:12000,:);
nCells = size(A,1);
%find all nearest neighnors and create sparse matrix
nnDist = 50;
[IDX, D] = rangesearch(A, A,nnDist);
IDX = cellfun(@(x,y) x(y>0), IDX,D,'uniformoutput',false);
pointsToMatch = arrayfun(@(x,y) repmat(x,y,1),1:numel(IDX), cellfun(@numel,IDX)','UniformOutput', false);
nnMatrix = sparse(cat(1,pointsToMatch{:}), cat(2,IDX{:})',1, nCells,nCells);

%%
load('/bigstore/GeneralStorage/Alon/Figures/R01June2019/DataForPlots.mat')
nReactions = 3;
nSpecies = 5;
%% Init cells, small frac infected
fracInf = 0.0001;
infected = rand(nCells,1)<fracInf;

fracMacs = 0.015;
Macs = rand(nCells,1)<fracMacs;

x= [~infected.*~Macs, infected.*~Macs, zeros(nCells,1), zeros(nCells,1), Macs]';
x = x(:);


nReactions = 3;
nSpecies = 5;

%
%x = DataForPlots.Macs05.xInit
% Plot Initial conditions
Macs = logical(x(5:nSpecies:end))
infected = logical(x(2:nSpecies:end))


MacCents = A(Macs,:);
searchRadiusVirus = 200;
searchRadiusMacs = 400;



[ids , dists] = rangesearch(A(logical(infected.*~Macs),:),MacCents,searchRadiusVirus);
MacsThatProduce = ~cellfun(@isempty, ids);

TNFofR=@(x) YukavaFit([1000, 1/150, 0], x);
[ids , dists] = rangesearch(MacCents(MacsThatProduce,:),A,searchRadiusMacs);

IC50=1;
amp = 24;
basln = 1;
TNFs = cellfun(@(x) sum(TNFofR(x)),dists);
TNFs(TNFs==inf)=0;
Lifetimes = (amp-basln)./(1+TNFs./IC50)+basln;


figure('color','w')
ax1 = axes( 'Units', 'Inches','position',[0.5 0.5 4.8 4.8]);
h(1) = scatter(A(logical(x(1:nSpecies:end)),2),A(logical(x(1:nSpecies:end)),3),30,1./Lifetimes(logical(x(1:nSpecies:end))),'filled','MarkerFaceAlpha', 0.2);
set(ax1,'xcolor','none','ycolor','none','CLim',[0,1]);
colormap('parula')
hold on;
h(2) = scatter(A(logical(x(2:nSpecies:end)),2),A(logical(x(2:nSpecies:end)),3),30,'r','filled','MarkerFaceAlpha', 0.8);
h(3) = scatter(A(logical(x(3:nSpecies:end)),2),A(logical(x(3:nSpecies:end)),3),30,'k','filled','MarkerFaceAlpha', 0.4);
h(4) = scatter(A(logical(x(4:nSpecies:end)),2),A(logical(x(4:nSpecies:end)),3),30,'g','filled','MarkerFaceAlpha', 0.6);
h(5) = scatter(A(logical(x(5:nSpecies:end)),2),A(logical(x(5:nSpecies:end)),3),30,'y','filled','MarkerFaceAlpha', 0.6);
%axis equal;

[hl, hlobj] = legend('Healthy','Infected','False positive','Dead post infection','Mac');
hl.Position = [0.8, 0.7, 0.14, 0.2];
hl.Box = 'off';
hl.FontSize = 14;
for i=6:10;
    hlobj(i).Children.MarkerSize=10;
end
cb = colorbar('Position',[0.8,0.15, 0.05, 0.4])
cb.Label.String = 'TNF\alpha';
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'Low','Medium','High'}
%%
figure(gcf)
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r300',['/bigstore/GeneralStorage/Alon/Figures/FiguresForLM/ModelWithMacs/' 'Macs0_5_0h']);
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/Figures/FiguresForLM/ModelWithMacs/' 'Macs0_5_h']);


%% Run model
xinitial=x;
VGR = 1/6;%viral growth rate
VI = 1/2; %viral infectivity

%savpath = ['/bigstore/GeneralStorage/Alon/Figures/FiguresForLM/ModelWithMacs/Macs0_5p_NewCornea_2']

[xOut, infDeathRate] = SSA_forSpatialSIDWithFixedMacs(x,nnMatrix,A,VGR,VI);

sum(reshape(xOut,nSpecies,[]),2)./sum(sum(reshape(xOut,nSpecies,[]),2))

%%
condName = 'Macs15_2'
DataForPlots.(condName).initDeathRate = 1./Lifetimes;
DataForPlots.(condName).xOut = xOut;
DataForPlots.(condName).xInit = xinitial;
DataForPlots.(condName).infDeathRate = infDeathRate;
DataForPlots.(condName).fracInf = fracInf;
DataForPlots.(condName).fracMacs = fracMacs;

%%
condName = 'Macs15_2'
x = DataForPlots.(condName).xOut;
infDeathRate = DataForPlots.(condName).infDeathRate;
figure('color','w')
ax1 = axes( 'Units', 'Inches','position',[0.5 0.5 4.8 4.8]);
h(1) = scatter(A(logical(x(1:nSpecies:end)),2),A(logical(x(1:nSpecies:end)),3),30,infDeathRate(logical(x(1:nSpecies:end))),'filled','MarkerFaceAlpha', 0.2);
set(ax1,'xcolor','none','ycolor','none','CLim',[0,1]);
colormap('parula')
hold on;
h(2) = scatter(A(logical(x(2:nSpecies:end)),2),A(logical(x(2:nSpecies:end)),3),30,'r','filled','MarkerFaceAlpha', 0.8);
h(3) = scatter(A(logical(x(3:nSpecies:end)),2),A(logical(x(3:nSpecies:end)),3),30,'k','filled','MarkerFaceAlpha', 0.4);
h(4) = scatter(A(logical(x(4:nSpecies:end)),2),A(logical(x(4:nSpecies:end)),3),30,'g','filled','MarkerFaceAlpha', 0.6);
h(5) = scatter(A(logical(x(5:nSpecies:end)),2),A(logical(x(5:nSpecies:end)),3),30,'y','filled','MarkerFaceAlpha', 0.6);
axis equal;

[hl, hlobj] = legend('Healthy','Infected','False positive','Dead post infection','Mac');
hl.Position = [0.8, 0.7, 0.14, 0.2];
hl.Box = 'off';
hl.FontSize = 14;
for i=6:10;
    hlobj(i).Children.MarkerSize=10;
end
            tl = title(['T = ' sprintf('%.1f',48) 'h']);
cb = colorbar('Position',[0.8,0.15, 0.05, 0.4])
cb.Label.String = 'TNF\alpha';
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'Low','Medium','High'}
%%
figure(gcf)
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r300',['/bigstore/GeneralStorage/Alon/Figures/FiguresForLM/ModelWithMacs/' 'Macs0_5_48h']);
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/Figures/FiguresForLM/ModelWithMacs/' 'Macs0_5_48h']);


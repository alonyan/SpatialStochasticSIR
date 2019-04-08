%% Make triangular grid
A = triangleGrid([-10 -10 10 10], [0,0], 1);
nCells = size(A,1);
%find all nearest neighnors and create sparse matrix
[IDX, D] = rangesearch(A, A,1.01);
IDX = cellfun(@(x,y) x(y>0), IDX,D,'uniformoutput',false);
pointsToMatch = arrayfun(@(x,y) repmat(x,y,1),1:numel(IDX), cellfun(@numel,IDX)','UniformOutput', false);
nnMatrix = sparse(cat(1,pointsToMatch{:}), cat(2,IDX{:})',1);

%% Init cells, small frac infected
fracInf = 0.005;
infected = rand(nCells,1)<fracInf;

x= [~infected, infected, zeros(nCells,1), zeros(nCells,1)]';
x = x(:);

figure('color','w')
ax1 = axes('position',[0,0,1,1])
scatter(A(logical(x(1:4:end)),1),A(logical(x(1:4:end)),2),[],'b','filled')
hold on;
scatter(A(logical(x(2:4:end)),1),A(logical(x(2:4:end)),2),[],'r','filled')
scatter(A(logical(x(3:4:end)),1),A(logical(x(3:4:end)),2),[],'k','filled')
scatter(A(logical(x(4:4:end)),1),A(logical(x(4:4:end)),2),[],'g','filled')
axis equal
set(ax1,'xcolor','w','ycolor','w')
shg
hl = legend('Healthy','Infected','False positive','Dead post infection')
hl.Location = 'northeastoutside';
hl.Box = 'off';

figure(gcf)
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-dpng','-r300',['/Users/Alonyan/Desktop/' 'Initial']);
%% Run model

VGR = 5000;%viral growth rate
VI = 1.5; %viral infectivity
alpha=10;%effect of TNF


xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,alpha)


%
figure('color','w')
ax1 = axes('position',[0,0,1,1])
scatter(A(logical(xOut(1:4:end)),1),A(logical(xOut(1:4:end)),2),[],'b','filled')
hold on;
scatter(A(logical(xOut(2:4:end)),1),A(logical(xOut(2:4:end)),2),[],'r','filled')
scatter(A(logical(xOut(3:4:end)),1),A(logical(xOut(3:4:end)),2),[],'k','filled')
scatter(A(logical(xOut(4:4:end)),1),A(logical(xOut(4:4:end)),2),[],'g','filled')
axis equal
set(ax1,'xcolor','w','ycolor','w')
an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['VI=' num2str(VI)], ['\alpha=' num2str(alpha)], ['VGR=' num2str(VGR)]},'LineStyle','none');

shg
hl = legend('Healthy','Infected','False positive','Dead post infection')
hl.Location = 'northeastoutside';
hl.Box = 'off';
%%
figure(gcf)
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-dpng','-r300',['/Users/Alonyan/Desktop/' 'HighAlpha2']);
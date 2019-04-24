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

% Run model

VGR = 1/6;%viral growth rate
VI = 1/2; %viral infectivity

%basalDeathRate=1/1000;%effect of TNF on deaths
%infDeathRate = 1/24;

basalDeathRate=1/100;%effect of TNF on deaths
infDeathRate = 1/4;

xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

pStats = sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2));

function pStats = modelRunForStats(nnMatrix,fracInf,VGR,VI,basalDeathRate,infDeathRate)
%% Init cells, small frac infected
nCells = size(nnMatrix,1);
infected = zeros(nCells,1);%<fracInf;
infected(116)=1; %always start with center pixel infected
x= [~infected, infected, zeros(nCells,1), zeros(nCells,1)]';
x = x(:);

% Run model

xOut = SSA_forSpatialSIDGrids(x,nnMatrix,[],VGR,VI,basalDeathRate,infDeathRate);

pStats = sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2));

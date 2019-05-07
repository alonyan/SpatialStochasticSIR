function [xOut, infDeathRate] = SSA_forSpatialSIDWithFixedMacs(x,neighMat,Grid,VGR,VI, varargin)
% Define propensity functions in terms of parameters, k, and states, x.
% each cell can be in 4 states - SIDF
w_k_x = @(x,k) [k(1)*x(1), k(2)*x(1), k(3)*x(2)];
% Specify stoichiometry
S = [ -1 -1  0 ;...
       1  0 -1 ;...
       0  1  0 ;
       0  0  1;
       0  0  0];%Macs don't play

%reactions: infection, false positive death, infected death

%number of cells
global nReactions nSpecies
nReactions = 3;
nSpecies = 5;

nC = numel(x)/nSpecies;

%make full stoichiometry matrix, block diagonal of single cell stoch mat
%SCell = repmat({S}, 1, nC);
%S = sparse(blkdiag(SCell{:}));

tstop = 48; %final time.

[xOut, infDeathRate] = Run_SSA(w_k_x,S,x,tstop,nC, neighMat,Grid,VGR,VI,varargin{:}); % call code to run stochastic simulation.

function [x, infDeathRate] = Run_SSA(prop_fun,S,x0,tstop,nC,neighMat,Grid,VGR,VI, varargin)


global nReactions nSpecies

t=0;
x = x0;     %% Specify initial conditions

arg.plot=ParseInputs('plot', false,varargin);
arg.save=ParseInputs('save', false,varargin);
arg.savpath=ParseInputs('savpath', '',varargin);
if arg.save && isempty(arg.savpath)
    error('If you want to save you must specify a savpath argument')
end


%Add Macs as sources of TNF
Macs = logical(x(5:nSpecies:end));
MacCents = Grid(Macs,:);
searchRadiusVirus = 200;
searchRadiusMacs = 400;
TNFofR=@(x) YukavaFit([1000, 1/150, 0], x);
IC50=1;
amp = 24;
basln = 1;
infDeathRateofTNF=@(x) 1./((amp-basln)./(1+x./IC50)+basln);
%scaling of infected/bystander death
p = [0.025 -0.0008];

%Init death rates
%find all macs that are close to infected cells, these would be producers
infected = logical(x(2:nSpecies:end));
[ids , ~] = rangesearch(Grid(logical(infected.*~Macs),:),MacCents,searchRadiusVirus);
MacsThatProduce = ~cellfun(@isempty, ids);
%calculate total TNF each cell sees
[~ , dists] = rangesearch(MacCents(MacsThatProduce,:),Grid,searchRadiusMacs);
TNFs = zeros(nC,1);
for ind1=find(~cellfun(@isempty, dists))'
    TNFs(ind1) = sum(TNFofR(dists{ind1}));
end


%calculate death rates due to this level of TNF
infDeathRate = infDeathRateofTNF(TNFs);
basalDeathRate=polyval(p,infDeathRate);%effect of TNF on deaths


%vector of virus infection. 0 for healthy or dead cells. initial infection=1 for infected cells.
v=x(2:nSpecies:end);

%Init params
k=zeros(1,nC*nReactions);
NNvLoad = neighMat*(v(:).*x(2:nSpecies:end));
k(1:nReactions:end) = VI./(1+1./NNvLoad);  %infection \propto virus load of NN
k(2:nReactions:end) = basalDeathRate;%basal death rate
k(3:nReactions:end) = infDeathRate;%./(1+1./v); %virus induced death \propto viral load of self

w=zeros(1,nC*nReactions);
for indCell=1:nC
    w((indCell-1)*nReactions+1:indCell*nReactions) = prop_fun(x(((indCell-1)*nSpecies+1):(indCell*nSpecies)), k(((indCell-1)*nReactions+1):(indCell*nReactions)));
end
%k(2:nReactions:end)=basalDeathRate;%basal death rate


if arg.plot
    figure('color','w')
    ax1 = axes( 'Units', 'Inches','position',[0.5 0.5 4.8 4.8]);
end

if arg.save
    outputVideo = VideoWriter(sprintf('%s%s',arg.savpath,'.avi'));
    outputVideo.FrameRate = 16;
    open(outputVideo);
end
frameNum=-1;

dists = {};
cellsInPlay = [];
while t<tstop

    if arg.plot
        if round(2*t)>frameNum
            cla
            frameNum=round(2*t);
            scatter(ax1,Grid(logical(x(1:nSpecies:end)),1),Grid(logical(x(1:nSpecies:end)),3),30,infDeathRate(logical(x(1:nSpecies:end))),'filled','MarkerFaceAlpha', 0.2);
            set(ax1,'xcolor','none','ycolor','none','CLim',[0,1]);
            colormap('parula')
            hold on;
            scatter(ax1,Grid(logical(x(2:nSpecies:end)),1),Grid(logical(x(2:nSpecies:end)),3),30,'r','filled','MarkerFaceAlpha', 0.8);
            scatter(ax1,Grid(logical(x(3:nSpecies:end)),1),Grid(logical(x(3:nSpecies:end)),3),30,'k','filled','MarkerFaceAlpha', 0.4);
            scatter(ax1,Grid(logical(x(4:nSpecies:end)),1),Grid(logical(x(4:nSpecies:end)),3),30,'g','filled','MarkerFaceAlpha', 0.6);
            scatter(ax1,Grid(logical(x(5:nSpecies:end)),1),Grid(logical(x(5:nSpecies:end)),3),30,'y','filled','MarkerFaceAlpha', 0.6);
            axis equal;
            tl = title(['T = ' sprintf('%.1f',frameNum/2) 'h']);
            %tl.Position(2)=10.5;
            
            if t==0
                [hl, hlobj] = legend('Healthy','Infected','False positive','Dead post infection','Mac');
                hl.Position = [0.8, 0.7, 0.14, 0.2];
                hl.Box = 'off';
                hl.FontSize = 14;
                for i=6:10;
                    hlobj(i).Children.MarkerSize=10;
                end
                
            %    an = annotation('textbox',[0.7, 0.1, 0.14, 0.5],'String',{['VI=' num2str(VI,2)], ['BasalDeathRate=' num2str(basalDeathRate,2)], ['InfectedDeathRate=' num2str(infDeathRate,2)], ['VGR=' num2str(VGR,2)]},'LineStyle','none', 'Fontsize', 12);
            end
            
            drawnow
            shg
            if arg.save
                frame = getframe(gcf);
                im = frame2im(frame);
                writeVideo(outputVideo,im);
            end
        end
    else
        if round(2*t)>frameNum
            frameNum=round(2*t)
        end
    end    
    
    
    
    
    
    w0 = sum(w);                               % sum of the prop. functions
    dt = 1/w0*log(1/rand);          %when's the next reaction?
    t = t+dt;                       % update time of next reaction, exp dist
    if t<=tstop
        r2w0=rand*w0;               % generate second random number and multiply by prop. sum
        
        i=find(cumsum(w)>r2w0,1,'first');
        
        whichReact = mod(i,nReactions);
        whichReact(whichReact==0)=nReactions;
        
        whichCell = floor(i/nReactions);
        whichCell(whichReact==nReactions) = whichCell-1;
        
        x((nSpecies*whichCell+1):nSpecies*(whichCell+1)) = x((nSpecies*whichCell+1):nSpecies*(whichCell+1)) + S(:,whichReact);% update the configuration
        
        if any(x<0)
            1+1
        end
        
        %x = x+S(:,i);                                 % update the configuration
        v = v+VGR*dt*x(2:nSpecies:end);                      % increase viral load for infected live cells
        
        
        
        
        %update production and death rates only if there's a new infection
        if whichReact==1 || whichReact==3
            infected = logical(x(2:nSpecies:end));
            [ids , ~] = rangesearch(Grid(logical(infected.*~Macs),:),MacCents,searchRadiusVirus);
            MacsThatProduce = ~cellfun(@isempty, ids);
            %calculate total TNF each cell sees
            [~ , dists] = rangesearch(MacCents(MacsThatProduce,:),Grid,searchRadiusMacs);
            TNFs = zeros(nC,1);
            for ind1=find(~cellfun(@isempty, dists))'
                TNFs(ind1) = sum(TNFofR(dists{ind1}));
            end
            %calculate death rates due to this level of TNF
            infDeathRate = infDeathRateofTNF(TNFs);
            basalDeathRate=polyval(p,infDeathRate);%effect of TNF on deaths
            k(2:nReactions:end) = basalDeathRate;%basal death rate
            k(3:nReactions:end) = infDeathRate;%./(1+1./v); %virus induced death \propto viral load of self
        end
        
        %update k
        NNvLoad = neighMat*(v(:).*x(2:nSpecies:end));
        k(1:nReactions:end) = VI./(1+1./NNvLoad);  %infection \propto virus load of NN
               
        %Propensity functions: all the probabilities of all the reactions
        %w = reshape(cell2mat(arrayfun(@(y,z) prop_fun(y{:},z{:}), mat2cell(x(:),nSpecies*ones(1,nC)), mat2cell(k(:),nReactions*ones(1,nC)),'uniformoutput', false))',1,[]);%sorry for the oneliner
        %w=zeros(1,nC*nReactions);
       % [find(neighMat*(v(:)))', find(x(3:nSpecies:end))']
       cellsInPlay = unique([cellsInPlay find(neighMat*(v(:)))', find(x(3:nSpecies:end))', find(x(4:nSpecies:end))', find(~cellfun(@isempty, dists))']);
        for indCell=cellsInPlay
            w((indCell-1)*nReactions+1:indCell*nReactions) = prop_fun(x(((indCell-1)*nSpecies+1):(indCell*nSpecies)), k(((indCell-1)*nReactions+1):(indCell*nReactions)));
        end
        
          
        
    end
    
end
if arg.save
    close(outputVideo);
end




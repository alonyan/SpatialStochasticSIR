function xOut = SSA_forSpatialSIDGrids(x,neighMat,Grid,VGR,VI,basalDeathRate,infDeathRate, varargin)
% Define propensity functions in terms of parameters, k, and states, x.
% each cell can be in 4 states - SIDF
w_k_x = @(x,k) [k(1)*x(1), k(2)*x(1), k(3)*x(2)];
% Specify stoichiometry
S = [ -1 -1  0 ;...
    1  0 -1 ;...
    0  1  0 ;
    0  0  1];

%reactions: infection, false positive death, infected death

%number of cells
global nReactions nSpecies
nReactions = 3;
nSpecies = 4;

nC = numel(x)/nSpecies;

%make full stoichiometry matrix, block diagonal of single cell stoch mat
SCell = repmat({S}, 1, nC);
S = sparse(blkdiag(SCell{:}));

tstop = 3*24; %final time.

xOut = Run_SSA(w_k_x,S,x,tstop,nC, neighMat,Grid,VGR,VI,basalDeathRate, infDeathRate,varargin{:}); % call code to run stochastic simulation.

function x = Run_SSA(prop_fun,S,x0,tstop,nC,neighMat,Grid,VGR,VI,basalDeathRate, infDeathRate, varargin)


global nReactions nSpecies

t=0;
x = x0;     %% Specify initial conditions

arg.plot=ParseInputs('plot', false,varargin);
arg.save=ParseInputs('save', false,varargin);
arg.savpath=ParseInputs('savpath', '',varargin);
if arg.save && isempty(arg.savpath)
    error('If you want to save you must specify a savpath argument')
end

%Init params
k=zeros(1,nC*nReactions);
k(2:nReactions:end)=basalDeathRate;%basal death rate

%vector of virus infection. 0 for healthy or dead cells. initial infection=1 for infected cells.
v=x(2:nSpecies:end);

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
while t<tstop
    
    if arg.plot
        if round(2*t)>frameNum
            frameNum=round(2*t);
            scatter(ax1,Grid(logical(x(1:nSpecies:end)),1),Grid(logical(x(1:nSpecies:end)),2),150,'b','filled');
            set(ax1,'xcolor','none','ycolor','none');
            hold on;
            scatter(ax1,Grid(logical(x(2:nSpecies:end)),1),Grid(logical(x(2:nSpecies:end)),2),150,'r','filled');
            scatter(ax1,Grid(logical(x(3:nSpecies:end)),1),Grid(logical(x(3:nSpecies:end)),2),150,'k','filled');
            scatter(ax1,Grid(logical(x(4:nSpecies:end)),1),Grid(logical(x(4:nSpecies:end)),2),150,'g','filled');
            axis equal;
            tl = title(['T = ' sprintf('%.1f',frameNum/2) 'h']);
            tl.Position(2)=10.5;
            
            if t==0
                [hl, hlobj] = legend('Healthy','Infected','False positive','Dead post infection');
                hl.Position = [0.8, 0.7, 0.14, 0.2];
                hl.Box = 'off';
                hl.FontSize = 14;
                for i=5:8;
                    hlobj(i).Children.MarkerSize=10;
                end
                
                an = annotation('textbox',[0.7, 0.1, 0.14, 0.5],'String',{['VI=' num2str(VI,2)], ['BasalDeathRate=' num2str(basalDeathRate,2)], ['InfectedDeathRate=' num2str(infDeathRate,2)], ['VGR=' num2str(VGR,2)]},'LineStyle','none', 'Fontsize', 12);
            end
            
            drawnow
            shg
            if arg.save
                frame = getframe(gcf);
                im = frame2im(frame);
                writeVideo(outputVideo,im);
            end
        end
    end
    
    %update k
    NNvLoad = neighMat*(v(:).*x(2:nSpecies:end));
    k(1:nReactions:end) = VI./(1+1./NNvLoad);  %infection \propto virus load of NN
    k(3:nReactions:end) = infDeathRate;%./(1+1./v); %virus induced death \propto viral load of self
    
    %Propensity functions: all the probabilities of all the reactions
    w = reshape(cell2mat(arrayfun(@(y,z) prop_fun(y{:},z{:}), mat2cell(x(:),nSpecies*ones(1,nC)), mat2cell(k(:),nReactions*ones(1,nC)),'uniformoutput', false))',1,[]);%sorry for the oneliner
    w0 = sum(w);                               % sum of the prop. functions
    dt = 1/w0*log(1/rand);          %when's the next reaction?
    t = t+dt;                       % update time of next reaction, exp dist
    if t<=tstop
        r2w0=rand*w0;               % generate second random number and multiply by prop. sum
        i=1;                                          % initialize reaction counter
        while sum(w(1:i))<r2w0             % what's the next reaction? increment counter until sum(w(1:i)) exceeds r2w0
            i=i+1;
        end
        x = x+S(:,i);                                 % update the configuration
        v = v+VGR*dt*x(2:nSpecies:end);                      % increase viral load for infected live cells
    end
    
end
if arg.save
    close(outputVideo);
end




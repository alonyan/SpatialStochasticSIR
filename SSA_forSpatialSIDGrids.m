function xOut = SSA_forSpatialSIDGrids(x,neighMat,Grid,VGR,VI,alpha)
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
nC = numel(x)/4;

%make full stoichiometry matrix, block diagonal of single cell stoch mat
SCell = repmat({S}, 1, nC);
S = sparse(blkdiag(SCell{:}));

tstop = 0.2; %final time.

xOut = Run_SSA(w_k_x,S,x,tstop,nC, neighMat,Grid,VGR,VI,alpha); % call code to run stochastic simulation.

function x = Run_SSA(prop_fun,S,x0,tstop,nC,neighMat,Grid,VGR,VI,alpha)
t=0;
x = x0;     %% Specify initial conditions

arg.plot=false;
%Init params
k=zeros(1,nC*3);
k(2:3:end)=alpha;%basal death rate

%vector of virus infection. 0 for healthy or dead cells. initial infection=1 for infected cells.
v=x(2:4:end);

if arg.plot
    figure('color','w')
    ax1 = axes('position',[0.1 0.1 0.8 0.8]);
end

while t<tstop  
       
if arg.plot  
    scatter(ax1,Grid(logical(x(1:4:end)),1),Grid(logical(x(1:4:end)),2),[],'b','filled');
    title(['T = ' num2str(t) '[au]']);   
    hold on;
    axis equal;    
    scatter(ax1,Grid(logical(x(2:4:end)),1),Grid(logical(x(2:4:end)),2),[],'r','filled');
    scatter(ax1,Grid(logical(x(3:4:end)),1),Grid(logical(x(3:4:end)),2),[],'k','filled');
    scatter(ax1,Grid(logical(x(4:4:end)),1),Grid(logical(x(4:4:end)),2),[],'g','filled');
    shg
    
    if t==0
        hl = legend('Healthy','Infected','False positive','Dead post infection');
        hl.Location = 'northeastoutside';
        hl.Box = 'off';
        set(ax1,'xcolor','w','ycolor','w');
        
        an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['VI=' num2str(VI)], ['\alpha=' num2str(alpha)], ['VGR=' num2str(VGR)]},'LineStyle','none');
    end

    drawnow
end   
    
    %update k
    k(1:3:end) = VI*neighMat*(v(:).*x(2:4:end));  %infection \propto virus load of NN               
    k(3:3:end) = v*(1+alpha); %virus induced death \propto viral load of self
    
    %Propensity functions: all the probabilities of all the reactions
    w = reshape(cell2mat(arrayfun(@(y,z) prop_fun(y{:},z{:}), mat2cell(x(:),4*ones(1,nC)), mat2cell(k(:),3*ones(1,nC)),'uniformoutput', false))',1,[]);%sorry for the oneliner
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
        v = v+VGR*dt*x(2:4:end);                      % increase viral load for infected live cells         
    end
end



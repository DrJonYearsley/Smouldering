%% A cellular automaton SIR model to look at propagation of smouldering combustion in peat

% Model was published in Belcher et al (2010) 
% "Baseline intrinsic flammability of Earth’s ecosystems estimated from 
% paleoatmospheric oxygen over the past 350 million years", PNAS Vol. 107,
% No 52, p22448-22453

clear all

% Summary of the experimental data from different O2 atmospheres
% 20.9% O2
%    Spread rate from peat ~ 0.001 cm/s (0.08 cm/min),
%    Burn duration ~165 mins
% 19% O2
%    Spread rate from peat ~ 1 mm/min (0.1 cm/min)
%    Burn duration ~115 mins
% 18% O2
%    Spread rate from peat ~ 0.00078 cm/s (0.0465 cm/min)
%    Burn duration ~84 mins

%% Model parameters
% The width and length of the burn box (in pixels)
n=50; % width
m=50; % length

dt = 0.2; % Time step in mins
dx = 0.2; % Physical step size of 1 pixel in cm (burn box is 10x10cm)
% So the fastest velocity that can be simulated is 0.2/0.1 = 1 cm/min

filepref = 'fireox_results';  % Prefix for filename of results

% Divide simulations between multiple files to reduce file size and memory usage
nRep = 1;    % Number of replicates per file (if ==1 then data is not saved)
nRep2 = 1;   % Number of files

tMax = 2000;  % Maximum number of iterations in a simulation
displayOn = true; % If true display simulation in figure

% Simulate across a range of beta and mu
%beta_a = [0.005 0.05];  % Min, max of beta parameter
%mu_a = [0.01,0.2]; % Min, max of mu parameter

% Simulate a precise value of beta and mu
beta_a = 0.0207 * [1 1]; % Min, max of beta parameter
mu_a = 0.0485 * [1 1]; % Min, max of mu parameter

%% Define the nearest neighbours of each cell

ind = (1:n*m)';
col = floor((ind-1)/n);
row = 1 + mod(ind-1,n);
% % Moore Nearest neighbour indicies
% nn = [n*col+ 1+mod(row-2,n), n*col+1+mod(row,n), ...
%     1+mod(ind-n-1,n*m),1+mod(ind+n-1,n*m)];

% % Von Neumann Nearest neighbour indicies
% nn = [n*col+1+mod(row-2,n), n*col+1+mod(row,n), ...
%     1+mod(ind-n-1,n*m),1+mod(ind+n-1,n*m), ...
%     1+mod(n*col+1+mod(row-2,n)-n-1,n*m),1+mod(n*col+1+mod(row-2,n)+n-1,n*m),...
%     1+mod(n*col+1+mod(row,n)-n-1,n*m), 1+mod(n*col+1+mod(row,n)+n-1,n*m)];

% von Neumann with no circular boundaries
% Last index is a dummy holding index
nn = (n*m+1)*ones(n*m+1,8);

dum = (row~=1);
nn([dum; false],1) = n*col(dum)+row(dum)-1;
dum = (row~=n);
nn([dum; false],2) = n*col(dum)+row(dum)+1;
dum = (col~=0);
nn([dum; false],3) = n*col(dum)+row(dum)-n;
dum = (col~=(m-1));
nn([dum; false],4) = n*col(dum)+row(dum)+n;
dum = (row~=1 & col~=0);
nn([dum; false],5) = n*col(dum)+row(dum)-n-1;
dum = (row~=n & col~=0);
nn([dum; false],6) = n*col(dum)+row(dum)-n+1;
dum = (row~=1 & col~=(m-1));
nn([dum; false],7) = n*col(dum)+row(dum)+n-1;
dum = (row~=n & col~=(m-1));
nn([dum; false],8) = n*col(dum)+row(dum)+n+1;

%%

log_beta_a = log(beta_a);
log_mu_a = log(mu_a);

for r2=1:nRep2
    if nRep>1  %%
        fileNum = 1;
        filepref_a = [filepref '_n' num2str(n) '_m' num2str(m)];
        filenameOut = [filepref_a '_#' num2str(fileNum) '.mat'];
        while exist(filenameOut,'file'),
            fileNum = fileNum+1;
            filenameOut = [filepref_a '_#' num2str(fileNum) '.mat'];
        end
    end
    
    burnDuration = zeros(nRep,1);
    burnProp = zeros(nRep,1);
    igniteTime = zeros(nRep,n*m);
    rowVel = zeros(nRep,n);
    rowDist = zeros(nRep,n);
    
    betaVal = exp(min(log_beta_a) + range(log_beta_a)*rand(nRep,1));
    muVal = exp(min(log_mu_a) + range(log_mu_a)*rand(nRep,1));
    
    for r = 1:nRep
        beta = betaVal(r);
        mu = muVal(r);
        
        if displayOn
            disp(['Beta = ' num2str(beta) ' and mu = ' num2str(mu)])
        end
        
        x = zeros(n,m);
        x(:,1) = 1; % Start simulation with a burning row
        x = [reshape(x,n*m,1); -1];
        
        tBurn = nan(n*m+1,1);
        tBurn(x==1) = 0;  % The initial burning cells
        tExtinct = nan(n*m+1,1);
        
        if nRep<2 && displayOn  % Display a pretty picture of the simulation
            hf1 = figure(1);
            set(hf1,'color',[1 1 1])
            
            hf2 = imagesc(1:n,1:m,reshape(x(1:(end-1)),n,m));
            set(gca,'xTickLabel',[],'YTickLabel',[],'XGrid','off','YGrid','off')
            
            caxis([-1 1])
            map = [0 0 0; 0.7895    0.4934    0.3142; 1 1 0];
            colormap(map)
            axis square
            
            disp('Hit any key to continue')
            pause
        end
        
        t = 1;
        keepGoing = true;
        while keepGoing && t<tMax
            toburn = (x==0);
            I = toburn & sum(beta*(x(nn)==1)>rand(n*m+1,8),2);
            R = (x==1) & mu>rand(n*m+1,1);
            
            tBurn(I) = t;
            tExtinct(R) = t;
            x(I) = 1;
            x(R) = -1;
            
            if nRep<2 && displayOn % Update the pretty picture of the simulation
                set(hf2,'CData',reshape(x(1:(end-1)),n,m))
                drawnow
                pause(0.01)
           end
            
            keepGoing = nnz(x==0)>0 && nnz(x==1)>0;
            t = t+1;
        end
        
        % If some material burning, add on the time until it finishes
        % burning, drawn from an exponential distribution
        burning = (x==1);
        
        tExtinct(burning) = t-log(rand(nnz(burning),1))/mu;
        
        keep = (1:n*m);
        burnDuration(r) = max(tExtinct(keep));
        burnt = x(keep)~=0;
        burnProp(r) = nnz(burnt)/n/m;
        for i=1:n
            rowDist(i) = max(col(burnt & row==i));
            rowVel(r,i) = rowDist(i) / tBurn(i+n*rowDist(i));
        end
        ind = ~isnan(tBurn);
        igniteTime(r,ind) = tBurn(ind)';
        
        if mod(r,500)==0,
            save(filenameOut,'mu','beta','n','m','nRep','tMax','betaVal','muVal',...
                'burnDuration','burnProp','rowVel','igniteTime','dt','dx','col','row')
        end
    end
    
    % Set nan velocities equal to zero
    ind = isnan(rowVel);
    rowVel(ind) = 0;
    
    igniteTime = sparse(igniteTime);
    rowVel = sparse(rowVel);
    
    if nRep>1
        % Variables for a line burn model
        save(filenameOut,'mu','beta','n','m','nRep','tMax','betaVal','muVal',...
            'burnDuration','burnProp','rowVel','igniteTime','dt','dx','col','row')
    end
end
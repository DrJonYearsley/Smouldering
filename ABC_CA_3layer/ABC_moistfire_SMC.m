% ABC_moistfire_SMC.m
%
% Script to perfrom ABC - Sequential Monte-Carlo
%
% This script applies seqABC to estimate parameter values for the
% model of smouldering in wet peat (3-layer cellular automata)
%
% This script calls  ABC_moistfire_simulation.m and ABC_distance.m.
% The script ABC_moistfire_setup.m is also required 
% 
% Sequential Monte-Carlo Algorithm
%  1. Set tolerance criterion to accept a simulation
%  2. Simulate picking parameters from a prior until tolerance reached.
%  3. Repeat steps 1 & 2 until n simulations (chains) created
%  4. Reduce tolerance (advance one time step)
%  5. Use simulations to create a prior
%  6. Return to step 1
% 
% This algorithm is described in Toni et al (2009), Approximate
% Bayesian computation scheme for parameter inference and model
% selection in dynamical systems, Interface, Vol. 6, p187-202.
% doi:10.1098/rsif.2008.0172
%
% Thus script uses a uniform prior distribution and a symetric proposal distribution so
% that acceptance along the MCMC is always 1.

%
% Written Jon Yearsley (jon.yearsley@ucd.ie) July 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

filepref = 'FIRE_ABC'; % Prefix to use for all output files
display = true;        % Print output of each Monte-Carlo iteration

chain.num = 5;    % Number of MCMC chains (particles) to run
chain.tsteps = 4; % Number of sequential steps to make
chain.epsilon = linspace(3,0.5,chain.tsteps); % Desired accuracy for each sequential step

% Names of variables to manipulate
states.names = {'sigmaI',  'alpha_bath', 'Hi','f_k', 'f_c', 'oxy_Thalf', 'oxy_A'};

% Define max an min values of the states. Prior is a uniform distribution
states.max = [5e-4, 5e-4,   2e4, 1.0, 1.0, 600, 100];
states.min = [5e-5, 5e-5, 1.5e4, 0.7, 0.7, 250, 5];

% Uniform prior distribution
states.prior = @(smax, smin) smin + (smax-smin).*rand(size(smax));
% Vector to save accept simulations
states.accept = nan(chain.tsteps, chain.num, length(states.max));
% Function to perform perturbation
states.perturb_sd = (states.max-states.min)/20;

% Gaussian perturbation kernel (perturb = find new parameter, K = transition probability kernel)
states.perturb = @(smean, ssd) smean+ssd.*randn(size(smean));
states.K = @(state1, smean, ssd, smax, smin) prod((sqrt(1/(2*pi))./ssd).*exp(-(state1-smean).^2./2/ssd.^2)./(abs(erf((smax-smean)./ssd/sqrt(2))/2)+abs(erf((smean-smin)./ssd/sqrt(2)))/2));

% % Uniform perturbation kernel (perturb = find new parameter, K = transition probability kernel)
% states.perturb = @(smean, ssd) smean + 2*ssd.*(rand(size(smean))-0.5); % Uniform perturbation kernel
% states.K = @(state1, smean, ssd, smax, smin) 1/prod((min(smax,smean+ssd)-max(smin,smean-ssd)));

range = (states.max-states.min)/10; % Default range

%% Set up default parameters for moistfire
P = ABC_moistfire_setup();
P.displayOn = false; % If true display simulation in figure
P.movieOn = false; % If true then save an animated gif movie of the simulation

% Define target values (spread rate cm/min & burn duration, mins)
% Data for 0.05% VMC
target.point = [0.06, 400]; % velocity [cm/min], tBurn [mins], temp [deg C]
target.tol = [0.01, 100]; % Velocity, tBurn, temp (+/- 1SD)
%  100% MC is vel=2.8cm/h, time=445mins

P.VMC_0 = 0.05; % Set moisture content
%% Find an unused filename    
fileNum = 1;
filepref_a = [filepref '_nx' num2str(P.n(1)) '_ny' num2str(P.n(2))];
filenameOut = [filepref_a '_#' num2str(fileNum) '.mat'];
while exist(filenameOut,'file'),
    fileNum = fileNum+1;
    filenameOut = [filepref_a '_#' num2str(fileNum) '.mat'];
end



%% Do sequential MCMC

chain.tries = zeros(chain.tsteps, chain.num); % Record the number of tires
chain.weights = ones(chain.tsteps, chain.num)/chain.num; % Record the weights
K = zeros(1,chain.num); % Vector to be filled with transition probs

sim = zeros(1,2);
out = nan(chain.num, 2);

for t=1:length(chain.epsilon)
    if display, 
        disp(['======== t=',num2str(t)]) 
    end
    
    for p=1:chain.num
        if display
            disp(['  ++++++++++ chain num=',num2str(p)])
        end
    
        dist = chain.epsilon(t)+1;
        while max(dist)>chain.epsilon(t) % Continue until a simulation is close to target (until max(dist)<=chain.epsilon)
            chain.tries(t,p) = chain.tries(t,p)+1;
            
            % Define initial states
            if t==1
                % Select particle parameters from prior distribution
                states.init = states.prior(states.max, states.min);
                weight_ind = 0;
            else
                % Select a particle proportional to the weights and generate a new starting value
                weight_ind = find((rand-eps)<=cumsum(chain.weights(t-1,:)),1);
                states.init = states.perturb(squeeze(states.accept(t-1, weight_ind, :))', states.perturb_sd);
                while any(states.init>states.max | states.init<states.min) % Make sure parameter stays within bounds
                    states.init = states.perturb(squeeze(states.accept(t-1,weight_ind,:))', states.perturb_sd);
                end
            end
            
            
            % Set the initial values for the parameters to be fitted
            for i=1:length(states.names)
               eval(['P.' states.names{i} '=states.init(i);']) 
            end
            
            % Perform simulation
            [times, rec] = ABC_moistfire_simulation(P);
            
            % Record results and compare with target values
            tburn=reshape(times.burn,P.n(1),P.n(2))*P.dt;
            ruler.n1 = [1:P.n(1)]*P.dx(1);
            ruler.n2 = [1:P.n(2)]*P.dx(2);
            
            n1 = repmat(ruler.n1',1,P.n(2));
            n2 = repmat(ruler.n2,P.n(1),1);
            n1 = n1(:);
            n2 = n2(:);
            tburn=tburn(:);
            ind = ~isnan(tburn);
            tmp=[ones(sum(ind),1) n1(ind) n2(ind)] \ tburn(ind);
            velx=1./tmp(2);
                       
            tmp = rec.T(:,3:end-3)+P.abs_zero;
            
            % If velocity and burn times are within sensible limits
            if velx<min(P.dx)/P.dt && max(times.extinct)<P.tMax && max(tmp(:))>100
                sim(1) = velx; % velocity
                sim(2) = myQuantile(tmp(:),0.95); % temp
                
                % Calculate distance bewteen taregt and simulation
                dist = ABC_distance(sim,target);
                
                if display
                    disp(sim)
                    disp([chain.tries(t,p) weight_ind dist' chain.epsilon(t)])
                end
            else
                dist = chain.epsilon(t)+1; % don't accept this simulation
            end
        end
        % Save accepted simulation outputs
        out(p,1) = sim(1); % velocity
        out(p,2) = sim(2); % temp
        
        
        list = strcat('P.', states.names);
        states.accept(t,p,:) = eval(['[ ' strjoin(list,' , ') ',]']); % Save new accepted state
        
        % Calculate the weight of the new particle
        K(p) = 0;
        if t>1
            for j=1:chain.num
                K(p) = K(p) + chain.weights(t-1,j) * states.K(squeeze(states.accept(t,p,:))', squeeze(states.accept(t-1,j,:))', states.perturb_sd, states.max, states.min);
            end
            % Normalize and save the new chain (particle) weights
            chain.weights(t,:) = 1./K; % Assumes a uniform prior
            chain.weights(t,:) = chain.weights(t,:) / sum(chain.weights(t,:)); % Normalise
        end
    end
    
    
    % Redefine perturbation to be proportional to variability between
    % particles
    states.perturb_sd = std( states.accept(t,:,:) )*2;
    
    ind = states.perturb_sd<1e-6;
    if nnz(ind)>0
        states.perturb_sd(ind) = range(ind);
    end
    
    save(filenameOut,'P','chain','states','out','target')
end


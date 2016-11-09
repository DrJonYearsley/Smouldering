% A reaction-diffusion cellular automaton model for smouldering peat
% where the cellular automata models heat diffusion
%
% Written: Jon Yearsley (jon.yearsley@ucd.ie) 9th July 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The cellular automaton has two layers:
%     1. A layer giving the temperature of each grid square
%     2. A layer giving the combustion state of each grid square
%
% All non-local effects are governed by heat diffusion in the temperature layer. 
% This is slightly different from standard CA which define a neighbourhood around
% each grid square. The dynamics of the temperature layer are driven by the heat equation with
% parameters suitable for peat. Boundary conditionas are that heat can be lost 
% through the sides and from the top of each grid square. The sides and top are 
% considered to be a heat bath of constant temperature. To try and be computationally efficient
% yet stable integration scheme the code uses operator splitting with an alternating direction implicit method
%
% There are 5 combustion states: unburned peat -> pyrolysis -> char -> oxidation -> ash
% The transition probabilities from one state to another are a function of the temperature 
% of the grid square. The function is sigmoidal with two parameters, Thalf and A
% The transition probabilities have the form
%    prob of transition = X / (1+X)
% where X = exp(A * (T-Thalf) / Thalf)
% T is the temperature of the grid square (Kelvin), 
% Thalf is the temperature when the transition probability is 0.5
% A is a parameter that determines how quickly the temperature causes the prob of transition
% to increase.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%close all

tic
% This is a development of the FIREOX model
% beta1 = pyrolysis due to heat from burning, I (S -> E) S being peat, E being char
% beta2 = pyrolysis due to heat from char, E (S -> E)
% alpha = onset of fire (E->I) I being burning char, 
% sigma = fire is extinguished (I->R), R being ash
% d = neighbourhood size of ignited peat

P.seed = 1; % if =0 then start with random seed, otherwise use value

P.n=[50, 50]; % width

P.dt = 0.1; % Time step in mins
P.dx = [0.5 0.5]; % Physical step size of 1 pixel in cm
% So the fastest velocity is 1 cell in 1 time step = 3/50/0.1 = 0.6 cm/min
% Shortest experiment duration ~ 30mins
% Peak mass loss reached at ~15mins for 35% O2. Max spread rate = 0.2cm/min
% tMax=20000 gives maximum duration of 2000mins


P.filepref = 'test';  % Prefix for filename of results
% Divide simulations between multiple files
% to reduce file size and memory usage
P.nRep = 1;    % Number of replicates per file (if ==1 then data is not saved)
P.nRep2 = 1;   % Number of files
P.tMax = 40000;  % Maximum number of iterations
P.displayOn = true; % If true display simulation in figure
P.movieOn = false; % If true then save an animated gif movie of the simulation
TmaxDisplay = 500;

% Fransden give 14.2MJ/kg (sd = 0.42) as heat emited by peat
% This is fairly independent of the moisture content
% Bomb calorimetry gives 19.4MJ/kg

P.Tb = 10; % [K] Temperature of surrounding heat bath
P.Ts = 50; % [K] Temperature of box sides (insulated)
P.Hi = 14.2*10^3; % Heat output of oxidising char [J/g] 
P.Hp = 0; % Heat output of pyrolysis [J/g]

P.Tinit = 600; % Temperature of igniter
P.tIgnition = 30/P.dt; % Time igniter is left on [time steps]

% Equation 2 from Benscoter et al (2011), 
% Vol 20, p 418 Int J of Wildland Fire
% Thermal conductivity of peat varies from 0.05*10^-2 to 0.15*10^-2 [W/cm/K]
P.k = 0.0015; % Thermal conductivity of peat [W/cm/K]
P.ks = 0; % Thermal conductivity of the sides [W/cm/K]
% 0.001 W/cm/K is conductivity of brick/wood

P.c = 0.92; % Specific heat of peat [J/g/K] (Oke 1997, Benscoter et al 2011)
% Thermal diffusivity = k/c/density [K/s /K cm^2 = cm^2/s]

% model parameters for transitions [Max and min values]
P.SC_Thalf_a = 180* [1 1];
P.SC_A_a = 50* [1 1];
P.CI_Thalf_a = 250* [1 1];
P.CI_A_a = 20* [1 1];
P.sigmaI_a = 0.03 * [1 1]; % Rate of oxidation per min
P.sigmaP_a = 0.02 * [1 1]; % Rate of pyrolysis per min
P.CV = 0.5; % If distribution allows variance to be set...

P.logScale = false; % If true select parameters from a log scale

% Define distribution of pyrolysis and oxidation times 
% Can be 'const','exp','unif'
stoch.pyro = 'norm'; 
stoch.oxy = 'norm'; 

% A factor governing heat losses throughout the sample to the heat bath
P.DAir = 0.002;  % Diffusivity constant from cell to heatbath


% Densities
% Peat = 151 +/- 3 kg /m^3 [10^-3 g/cm^3]
% Ash  = 36 +/- 2 kg /m^3  [10^-3 g/cm^3]
% Char = 189 +/- 4 kg /m^3 [10^-3 g/cm^3]
density.peat = 151 *10^-3; % g cm^-3
density.char = 189*10^-3; % g cm^-3
density.ash = 36*10^-3; % g cm^-3


state = struct('S',0,'P',1,'C',2,'I',3,'R',-1);

%% Set random number generator seed
if P.seed==0
    rng('shuffle')
    P.seed = rng;
else
    rng(P.seed)
    P.seed = rng;
end
%% Split operator

alpha = (P.k/P.c/density.peat * 60) * P.dt ./ P.dx.^2; % thermal diffusivity [L^2 / T]
alphas = (P.ks/P.c/density.peat * 60) * P.dt ./ P.dx.^2; % thermal diffusivity [L^2 / T]
alphaa = (P.DAir * 60) * P.dt; % thermal diffusivity [L^2 / T]


A1_setup = [
    1, 1, -(alpha(1)+alphas(1)+alphaa)/2;  % Diagonal elements
    (2:P.n(1)-1)', (2:P.n(1)-1)', repmat(-(2*alpha(1)+alphaa), P.n(1)-2,1)/2; 
    P.n(1), P.n(1), -(alpha(1)+alphas(1)+alphaa)/2;
    (1:P.n(1)-1)',(2:P.n(1))', repmat(alpha(1)/2, P.n(1)-1,1); % Off-diagonal elements
    (2:P.n(1))',(1:P.n(1)-1)', repmat(alpha(1)/2, P.n(1)-1,1); % Off-diagonal elements
    ];

A2_setup = [
    1, 1, -(alpha(2)+alphas(2)+alphaa)/2;  % Diagonal elements
    (2:P.n(2)-1)', (2:P.n(2)-1)', repmat(-(2*alpha(2)+alphaa), P.n(2)-2,1)/2; 
    P.n(2), P.n(2), -(alpha(2)+alphas(2)+alphaa)/2;
    (1:P.n(2)-1)',(2:P.n(2))', repmat(alpha(2)/2, P.n(2)-1,1); % Off-diagonal elements
    (2:P.n(2))',(1:P.n(2)-1)', repmat(alpha(2)/2, P.n(2)-1,1); % Off-diagonal elements
    ];

% Create explicit update matrices along dimension 1 and 2
A1_ex = sparse(A1_setup(:,1), A1_setup(:,2), A1_setup(:,3), P.n(1), P.n(1));
A2_ex = sparse(A2_setup(:,1), A2_setup(:,2), A2_setup(:,3), P.n(2), P.n(2));

% Create implicit update matrices along dimension 1 and 2
A1_imp = speye(P.n(1)) - A1_ex;
A2_imp = speye(P.n(2)) - A2_ex;

% Create boundary terms
b = P.Tb *alphaa * ones(P.n);
b(1,:) = b(1,:) + P.Ts*alphas(1);
b(end,:) = b(end,:) + P.Ts*alphas(1);
b(:,1) = b(:,1) + P.Ts*alphas(2);
b(:,end) = b(:,end) + P.Ts*alphas(2);
b = sparse(b)/2;

%%
% Define stochatic functions for pyrolysis and oxidation times
if strcmp(stoch.oxy,'exp')
  stoch.oxyFn = @(m,dim) -log(rand(dim))/m(1);
elseif strcmp(stoch.oxy,'unif')
  stoch.oxyFn = @(m,dim) 2*rand(dim)/m(1);
elseif strcmp(stoch.oxy,'gamma')
  stoch.oxyFn = @(m,dim) random(makedist('gamma','a',(m(2)/m(1))^2,'b',m(1)/m(2)^2),dim);
elseif strcmp(stoch.oxy,'norm')
  stoch.oxyFn = @(m,dim) max(5,1/m(1)+1/m(2)*randn(dim));
else
  stoch.oxyFn = @(m,dim) ones(dim)/m(1);
end

if strcmp(stoch.pyro,'exp')
  stoch.pyroFn = @(m,dim) -log(rand(dim))/m(1);
elseif strcmp(stoch.pyro,'unif')
  stoch.pyroFn = @(m,dim) 2*rand(dim)/m(1);
elseif strcmp(stoch.pyro,'gamma')
  stoch.pyroFn = @(m,dim) random(makedist('gamma','a',(m(2)/m(1))^2,'b',m(1)/m(2)^2),dim);
elseif strcmp(stoch.pyro,'norm')
  stoch.pyroFn = @(m,dim) max(5,1/m(1)+1/m(2)*randn(dim));
else
  stoch.pyroFn = @(m,dim) ones(dim)/m(1);
end
%%

log_sigmaI_a = log(P.sigmaI_a);
log_sigmaP_a = log(P.sigmaP_a);
log_SC_Thalf_a = log(P.SC_Thalf_a);
log_SC_A_a = log(P.SC_A_a);
log_CI_Thalf_a = log(P.CI_Thalf_a);
log_CI_A_a = log(P.CI_A_a);

for r2=1:P.nRep2
    if P.nRep>1  %%
        fileNum = 1;
        filepref_a = [P.filepref '_n' num2str(n) '_m' num2str(m)];
        filenameOut = [filepref_a '_#' num2str(fileNum) '.mat'];
        while exist(filenameOut,'file'),
            fileNum = fileNum+1;
            filenameOut = [filepref_a '_#' num2str(fileNum) '.mat'];
        end
        save(filenameOut,'P','density','stoch','col','row')

    end

    times = struct('pyro',zeros(prod(P.n),P.nRep),'char',zeros(prod(P.n),P.nRep), ...
        'burn',zeros(prod(P.n),P.nRep),'extinct',zeros(prod(P.n),P.nRep));
    
    if P.logScale
        sigmaIVal = exp(min(log_sigmaI_a) + (max(log_sigmaI_a)-min(log_sigmaI_a))*rand(P.nRep,1));
        sigmaPVal = exp(min(log_sigmaP_a) + (max(log_sigmaP_a)-min(log_sigmaP_a))*rand(P.nRep,1));
        SC.ThalfVal = exp(min(log_SC_Thalf_a) + (max(log_SC_Thalf_a)-min(log_SC_Thalf_a))*rand(P.nRep,1));
        SC.AVal = exp(min(log_SC_A_a) + (max(log_SC_A_a)-min(log_SC_A_a))*rand(P.nRep,1));
        CI.ThalfVal = exp(min(log_CI_Thalf_a) + (max(log_CI_Thalf_a)-min(log_CI_Thalf_a))*rand(P.nRep,1));
        CI.AVal = exp(min(log_CI_A_a) + (max(log_CI_A_a)-min(log_CI_A_a))*rand(P.nRep,1));
    else
        sigmaIVal = min(P.sigmaI_a) + (max(P.sigmaI_a)-min(P.sigmaI_a))*rand(P.nRep,1);
        sigmaPVal = min(P.sigmaP_a) + (max(P.sigmaP_a)-min(P.sigmaP_a))*rand(P.nRep,1);
        SC.ThalfVal = min(P.SC_Thalf_a) + (max(P.SC_Thalf_a)-min(P.SC_Thalf_a))*rand(P.nRep,1);
        SC.AVal = min(P.SC_A_a) + (max(P.SC_A_a)-min(P.SC_A_a))*rand(P.nRep,1);
        CI.ThalfVal = min(P.CI_Thalf_a) + (max(P.CI_Thalf_a)-min(P.CI_Thalf_a))*rand(P.nRep,1);
        CI.AVal = min(P.CI_A_a) + (max(P.CI_A_a)-min(P.CI_A_a))*rand(P.nRep,1);
    end
      
    
    for r = 1:P.nRep
        sigmaI = sigmaIVal(r);
        sigmaP = sigmaPVal(r);
        SC_Thalf = SC.ThalfVal(r);
        CI_Thalf = CI.ThalfVal(r);
        SC_A = SC.AVal(r);
        CI_A = CI.AVal(r);

        
        % Predefine some quantities
        ef.SC1 = exp(-SC_A);
        ef.SC2 = SC_A/SC_Thalf;
        ef.CI1 = exp(-CI_A);
        ef.CI2 = CI_A/CI_Thalf;
        expfun = zeros(prod(P.n),1);
 
        % Set the initial temperature
        T = P.Ts * ones(P.n);
        T(:,1) = P.Tinit;
               
        x = zeros(P.n);
        x(:,1) = state.I; % Start simulation with a burning row
        
        % Pre-assign times for pyrolysis and oxidation
        pyroTime = stoch.pyroFn((sigmaP*P.dt)*[1, 1/P.CV], P.n);
        oxyTime = stoch.oxyFn((sigmaI*P.dt)*[1, 1/P.CV], P.n);
        tPyro = nan(prod(P.n),1);  % Time at which cell starts pyrolysis
        tBurn = nan(prod(P.n),1);
        tBurn(x==state.I) = 0;  % The initial burning cells
        tExtinct = nan(prod(P.n),1);
        tExtinct(x==state.I) = stoch.oxyFn(sigmaI*P.dt*[1, 1/P.CV], [nnz(x==state.I),1]);
        tChar = nan(prod(P.n),1);
        data.temp = zeros(P.tMax, 10);
        step = [0:4]*P.n(1)*5;
        data.tempind = [15+step, 35+step];
        
        if P.nRep<2 && P.displayOn  % Display a pretty picture of the simulation
            hf1 = figure(1);
            set(hf1,'color',[1 1 1])
            
            hf2 = imagesc(1:P.n(1),1:P.n(2),x);
            set(gca,'xTickLabel',[],'YTickLabel',[],'XGrid','off','YGrid','off')
            t1= title('t = 0.0 mins','FontSize',15);
            caxis([-1 3])
            % Colours for ash, oxydation, char, pyrolysis and peat
            map = [0.98 0.98 0.98 ;  0.7895    0.4934    0.3142;  1 1 0 ;  0.2 0.2 0.2; 1 0.2 0.2];
            colormap(map)
            axis square

            hf3 = figure(2);
            set(hf3,'color',[1 1 1])
            
            hf4 = imagesc(1:P.n(1),1:P.n(2), T);
            set(gca,'xTickLabel',[],'YTickLabel',[],'XGrid','off','YGrid','off')
            axis square
            caxis([0 TmaxDisplay])
            t2= title('Temperature distribution, t = 0.0 mins','FontSize',15);
            colorbar
            
            disp('Hit any key to continue')
            pause
            if P.movieOn
                f = getframe(hf1);
                [im1,map1] = rgb2ind(f.cdata,256,'nodither');
                im1(1,1,1,P.tMax) = 0;
            end
        end
        
        t = 1;
        keepGoing = true;
        while keepGoing && t<P.tMax
            sup = (x==state.S);
            cha = (x==state.C);
            pyro = (x==state.P);
            oxy = (x==state.I);
            
            E = false(P.n);
            I = E;
            R = E;
            C = E;
 
            expfun(sup) = ef.SC1 * exp(ef.SC2*T(sup));
            E(sup) = expfun(sup)./(1+expfun(sup))>rand(nnz(sup),1);
            C(pyro) = t>tChar(pyro);            
            expfun(cha) = ef.CI1 * exp(ef.CI2*T(cha));
            I(cha) =  expfun(cha)./(1+expfun(cha))>rand(nnz(cha),1);
            R(oxy) = t>tExtinct(oxy);
           
            tPyro(E) = t;
            tChar(E) = t + pyroTime(E);
            tBurn(I) = t;
            tExtinct(I) = t + oxyTime(I);
            
            x(E) = state.P;
            x(C) = state.C;
            x(I) = state.I;
            x(R) = state.R;
                        
            % Calculate heat sources and sinks
            S = zeros(P.n);
            S(pyro) = P.Hp/P.c./(tChar(pyro)-tPyro(pyro));
            S(oxy) = P.Hi/P.c./(tExtinct(oxy)-tBurn(oxy));
            s = sparse(S)/2; % Divide by 2 for split operator calculations
            
            %% Heat transport
            % Update temp using implict update
            % Step 1
            %Tnew = T * A21; % Explicit along dimension 2
            T = A1_imp \ (T + (T * A2_ex) + b + s); % Implicit along dimension 1
            % Step 2
            %Tnew = A21 * T; % Explicit along dimension 1
            T = (T + (A1_ex * T) + b + s) / A2_imp; % Implicit along dimension 2
             
            data.temp(t,:) = T(data.tempind);
            if t<P.tIgnition
               T(:,1) = max([T(:,1),P.Tinit*ones(P.n(1),1)],[],2);
            end
            
            if P.nRep<2 && P.displayOn && floor(mod(t,1/P.dt))==0 % Update the pretty picture of the simulation
                set(t1,'String',['t=' num2str(round(t*P.dt)) 'mins'])
                set(t2,'String',['Temperature Distribution, t=' num2str(round(t*P.dt)) 'mins'])
                set(hf2,'CData', x)
                set(hf4,'CData', T)
                drawnow
                
                if P.movieOn
                    f = getframe(hf1);
                    im1(:,:,1,t) = rgb2ind(f.cdata,map1,'nodither');
                end
                
            end
            
            % Don't stop until there are no more burning cells
            if any(x(:)==state.S | x(:)==state.P | x(:)==state.C) && any(x(:)==state.I)
                keepGoing = true;
            else
                keepGoing=false;
            end
            t = t+1;
        end
        if P.displayOn && P.movieOn && t<P.tMax
            for t2 = t+1:tMax
                im1(:,:,1,t2) = rgb2ind(f.cdata,map1,'nodither');
            end
        end       
        keep = (1:prod(P.n));

        % Save the important times from the simulation
        times.pyro(:,r) = tPyro(keep);
        times.char(:,r) = tChar(keep);
        times.burn(:,r) = tBurn(keep);
        times.extinct(:,r) = tExtinct(keep);

        if mod(r,500)==0,
            save(filenameOut,'P','state','density','stoch','col','row', ...
                'sigmaI_a','sigmaP_a', 'SC_Thalf_a', 'SC_A_a', 'CA_Thalf_a', 'CI_A_a', ...
                'sigmaIVal','sigmaPVal','SC','CI','times')
         end
    end

    data.temp = data.temp(1:t,:);
    
    if P.nRep>1
        % Variables for a line burn model
        save(filenameOut,'P','state','density','stoch','col','row', ...
            'sigmaI_a','sigmaP_a', 'SC_Thalf_a', 'SC_A_a', 'CA_Thalf_a', 'CI_A_a', ...
            'sigmaIVal','sigmaPVal','SC','CI','times')
    end
end

if P.displayOn && P.movieOn
    % Save movie of the simulation
    imwrite(im1,map1,'fireox.gif','DelayTime',0,'LoopCount',0) %g443800
end
toc

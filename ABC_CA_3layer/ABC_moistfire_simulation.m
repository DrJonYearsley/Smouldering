function [times, rec] = ABC_moistfire_simulation(P)
% function [times, rec] = ABC_moistfire_simulation(P)
%
% A function that starts a simulation of smouldering in moist peat 
% using the parameters in structure P.
%
% This script calls ABC_moistfire_setup.m to initialize parameters
%
% Inputs:
%   P     A structure contaiing the parameters for the simulation
%
% Outputs:
%   times A structure containing the simulations times when each
%         grid square first entered a combustion state (possible
%         states are:
%            times.pyro=pyrolysis
%            times.char=char (end of pyrolysis)
%            times.burn=oxidation
%            times.extinct=ash (combustion has ended)
%
%   rec   A structure that contains the temperature and moisture
%         data from the simulation
%            rec.T(t,:)  temperature data for time t in the simulation
%            rec.v(t,:)  moisture data for time t in the simulation
%
% Written Jon Yearsley (jon.yearsley@ucd.ie) July 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up states and arrays for recording
% Different states: Fuel, Pyrolysis, Char, Oxidation, Ash
state = struct('F',0,'P',1,'C',2,'O',3,'A',-1);

% The temperatures from the simulation
rec.T = nan(floor(P.tMax/P.dt),length(P.indRec));  
% The moisture content from the simulation
rec.v = nan(floor(P.tMax/P.dt),length(P.indRec));

%% Solution to moisture dynamics
% Gravimetric moisture content, GMC
% Volumetric moisture content, VMC
%
% GMC = m_water / [m_water(t=0) + m_peat]
% GMC_dry = m_water / P.density_peat
% m_water = mass of water per kg of soil
%         = VMC * P.density_water
% P.density_peat = mass [kg] of dry peat per m^3 of soil
% m_peat = mass [kg] of dry peat per kg of soil
%
% VMC = GMC_dry * P.density_peat / P.density_water
%
% theta = m_water(t)/m_water(t=0)
%       = VMC(t) / VMC_0
%               m_water(t) = m(t) - m_peat = m(t)-m(t=infinity)
%               VMC_0 is the starting VMC
%
% d_theta/dt = -A exp(-E/(RT)) theta^n
% The solution is
% theta(Delta_t)^-(n-1) = theta(0)^-(n-1) + (n-1) * A * exp[-E/(R T)] Delta_t
% 1/VMC(Delta_t)^(n-1) = 1/VMC(0)^(n-1) + (n-1)*A*exp[-E/(R T)] Delta_t / VMC_0^(n-1)

% Parameter values from Chen et al 2011, Energy Fuels, Vol 25, p797-803
% E = 68500.0 J/mol
% A = exp(8.30) /sec
% n = 2.5
% R = 8.3144621 J/K/mol
%

% Model from Filkov et al (2011) Energy and Fuels
% phi = relative humidity
% B = multiplier = 0.46810^4  -  0.49*10^5 [K^0.5 s^-1]
% L = activation energy = 36029 - 42346 [J / mol]
% Model is
% dm/dt = - B (m-m_inf) exp(-L/R/T) * (1-phi) /sqrt(T)
% where m_inf = mass of sample at time infinity (i.e. mass of dry peat)
%       m = mass of sample at time t (peat and water)

if strcmp(P.drying_model,'chen');
    exp_dry = exp(-P.E_dry_chen/P.R/400);
    dry = @(T, Delta_t, v_0) (v_0.^(-3/2) + P.A_dry_chen * exp_dry.^(400./T).*Delta_t).^(-2/3);
else
    exp_dry=exp(-P.E_dry_filkov/P.R/400);
    dry = @(T, Delta_t, v_0) v_0.*exp(-P.A_dry_filkov*exp_dry.^(400./T).*Delta_t*(1-P.phi)./sqrt(T));
end

%% Define initial moisture contents, temperatures and states

ind = (1:prod(P.n))';
ind_x_lower = 1+(0:P.n(2)-1)*P.n(1);
ind_x_upper = P.n(1)+(0:P.n(2)-1)*P.n(1);
ind_y_lower = 1:P.n(1);
ind_y_upper = prod(P.n)-P.n(1)+(1:P.n(1));
ind_x_inner = ind;
ind_y_inner = ind;
ind_x_inner([ind_x_upper,ind_x_lower]) = [];
ind_y_inner([ind_y_upper, ind_y_lower]) = [];

trans = reshape(reshape(ind,P.n(1),P.n(2))',prod(P.n),1); % Use to transpose matrix
invtrans = reshape(reshape(ind,P.n(2),P.n(1))',prod(P.n),1); % Use to transpose matrix back

% Set the initial temperature
T = P.Tinit * ones(prod(P.n),1);
if P.IgniterConfig ==1
    ind_igniter = ind_x_lower;
elseif P.IgniterConfig ==2
    ind_igniter = floor(P.n(1)/2)+P.n(1)*([1:P.n(2)]-1);
elseif P.IgniterConfig==3
    ind_igniter = floor(P.n(2)/2)*P.n(1)+[1:P.n(1)];
elseif P.IgniterConfig==4
    ind_igniter = floor(P.n(2)/2)*P.n(1)+ceil(P.n(1)/2);
elseif P.IgniterConfig ==5
    ind_igniter = ind_y_lower;
elseif P.IgniterConfig ==6
    ind_igniter = ind_y_upper;
else
    ind_igniter = ind_x_upper;
end
T(ind_igniter) = P.Tigniter;

sponemat = speye(prod(P.n));


% Set the initial stat
x = zeros(prod(P.n),1);
x(ind_igniter) = state.O; % Start simulation with a burning row

% Set the initial volumetric moisture content
% Let v = VMC/VMC_0 so that initially all grid cells have values 1 (except
% the heat baths)
v = ones(size(T));
v_new = v;

% Set timestep to be half of P.dt (for operator splitting) and convert mins into seconds
% Factor of 60 converts from mins (P.dt units) into seconds
delta_t = 60*P.dt/2;


%% Split operator
% For the moment assume that the A operator acts as if T_new = A*T

% Setup matrices, column 1 gives first index, column 2 gives second index,
% column 3 gives the value of the matrix element

fac = 1./P.dx.^2;
% A1 is dimension 1 update matrix
A1_setup = zeros(3*prod(P.n),3); % this will be a tri-diagonal matrix
% Diagonal elements
A1_setup(ind,1:2) = [ind,ind];
% Off-diagonal upper elements
A1_setup(prod(P.n)+ind,1:2) = [ind,ind+1];
A1_setup(prod(P.n)+ind_x_upper,1) = 0;
% Off-diagonal lower elements
A1_setup(2*prod(P.n)+ind,1:2) = [ind,ind-1];
A1_setup(2*prod(P.n)+ind_x_lower,1) = 0;

% A2 is dimension2 update matrix (to keep this matrix tridiagonal (faster) is operates on the transpose of the T matrix)
A2_setup = zeros(3*prod(P.n),3); % this will be a tri-diagonal matrix (has to operate on transpose of T matrix)
% Diagonal elements
A2_setup(ind,1:2) = [ind,ind];
% For A2_setup we can replace ind_y_upper by trans(ind_x_upper) if region is square
% Off-diagonal upper elements
A2_setup(prod(P.n)+ind,1:2) = [ind,ind+1];
A2_setup(prod(P.n)+invtrans(ind_y_upper),1) = 0;
% Off-diagonal lower elements
A2_setup(2*prod(P.n)+ind,1:2) = [ind,ind-1];
A2_setup(2*prod(P.n)+invtrans(ind_y_lower),1) = 0;

% Remove extraneous elements from off-diagonals
id1 = A1_setup(:,1)>0;
id2 = A2_setup(:,1)>0;

%% Setup sparse matrices for A1 & A2 explicit and implicit, b1 and b2
if ~strcmp(computer,'MACI64')
    A1_ex_mask = spinit(A1_setup(id1,1), A1_setup(id1,2),  prod(P.n)*[1,1]);
    A2_ex_mask = spinit(A2_setup(id2,1), A2_setup(id2,2),  prod(P.n)*[1,1]);
    b1_mask = spinit(1:prod(P.n),ones(prod(P.n),1),[prod(P.n),1]);
end

%%
% Define stochatic functions for pyrolysis and oxidation times
if strcmp(P.oxyDist,'exp')
    oxyFn = @(m,dim) -log(rand(dim))/m(1);
elseif strcmp(P.oxyDist,'unif')
    oxyFn = @(m,dim) 2*rand(dim)/m(1);
elseif strcmp(P.oxyDist,'gamma')
    oxyFn = @(m,dim) random(makedist('gamma','a',(m(2)/m(1))^2,'b',m(1)/m(2)^2),dim);
elseif strcmp(P.oxyDist,'norm')
    oxyFn = @(m,dim) max(5,1/m(1)+1/m(2)*randn(dim));
else
    oxyFn = @(m,dim) ones(dim)/m(1);
end

if strcmp(P.pyroDist,'exp')
    pyroFn = @(m,dim) -log(rand(dim))/m(1);
elseif strcmp(P.pyroDist,'unif')
    pyroFn = @(m,dim) 2*rand(dim)/m(1);
elseif strcmp(P.pyroDist,'gamma')
    pyroFn = @(m,dim) random(makedist('gamma','a',(m(2)/m(1))^2,'b',m(1)/m(2)^2),dim);
elseif strcmp(P.pyroDist,'norm')
    pyroFn = @(m,dim) max(5,1/m(1)+1/m(2)*randn(dim));
else
    pyroFn = @(m,dim) ones(dim)/m(1);
end

t = 1;
keepGoing = true; 


%%

times = struct('pyro',zeros(prod(P.n),1),'char',zeros(prod(P.n),1), ...
    'burn',zeros(prod(P.n),1),'extinct',zeros(prod(P.n),1));

sigmaI = P.sigmaI;
sigmaP = P.sigmaP;
pyro_Thalf = P.pyro_Thalf;
oxy_Thalf = P.pyro_Thalf;
pyro_A = P.pyro_A;
oxy_A = P.oxy_A;


% Pre-assign times for pyrolysis and oxidation
pyroTime = pyroFn((sigmaP)*[1, 1/P.CV], [prod(P.n),1])*P.density_peat * prod(P.dx)/P.dt; % [time steps]
oxyTime = oxyFn((sigmaI)*[1, 1/P.CV], [prod(P.n),1])*P.density_peat * prod(P.dx)/P.dt;   % [time steps]
tPyro = nan(prod(P.n),1);  % Time at which cell starts pyrolysis
tBurn = nan(prod(P.n),1);
tBurn(x==state.O) = 0;  % The initial burning cells
tExtinct = nan(prod(P.n),1);
tExtinct(x==state.O) = oxyTime(x==state.O); % The initial burning cells
tChar = nan(prod(P.n),1);

% Heat source units J / s / cm3 of peat
heatsource_pyro = P.Hp * P.density_peat ./ pyroTime / P.dt/60;
heatsource_oxy = P.Hi * P.density_peat./ oxyTime / P.dt/60;

% Predefine some quantities
if pyro_A>0
    ef.SC2 = log(3)/pyro_A;
    ef.SC1 = 3*exp((pyro_Thalf)/pyro_A);
end
if oxy_A>0
    ef.CI2 = log(3)/oxy_A;
    ef.CI1 = 3*exp((oxy_Thalf)/oxy_A);
end
expfun = zeros(prod(P.n),1);

while keepGoing && t<=P.tMax/P.dt
    
    % Record temp and moisture
    rec.T(t,:) = T(P.indRec);
    rec.v(t,:) = P.VMC_0 * v(P.indRec);

    fuel = (x==state.F);
    cha = (x==state.C);
    pyro = (x==state.P);
    oxy = (x==state.O);
    
    % Keep track of states that have newly changed
    PY = false(prod(P.n),1);
    C = false(prod(P.n),1);
    O = false(prod(P.n),1);
    A = false(prod(P.n),1);
    
    if pyro_A>0
        expfun(fuel) = ef.SC1 * exp(-ef.SC2*T(fuel));
        PY(fuel) = expfun(fuel)./(1+expfun(fuel))<rand(nnz(fuel),1);
    else
        PY(fuel) = T(fuel)>pyro_Thalf;
    end
    C(pyro) = t>tChar(pyro);
    if oxy_A>0
        expfun(cha) = ef.CI1 * exp(-ef.CI2*T(cha));
        O(cha) =  expfun(cha)./(1+expfun(cha))<rand(nnz(cha),1);
    else
        O(cha) = T(cha)>oxy_Thalf;
    end
    A(oxy) = t>tExtinct(oxy);
    
    % Update states that are due to change
    tPyro(PY) = t;
    tChar(PY) = t + pyroTime(PY);
    tBurn(O) = t;
    tExtinct(O) = t + oxyTime(O);
    
    x(PY) = state.P;
    x(C) = state.C;
    x(O) = state.O;
    x(A) = state.A;
    
    %% Update moisture
    indEvap = T>P.Tmin_evap & v>P.VMC_min/P.VMC_0;
    v_new(indEvap) = dry(T(indEvap), 2*delta_t, v(indEvap));
    dvdt  = (v_new-v)/(2*delta_t);
    k = P.k_peat + P.k_water*P.f_k * (v+v_new)/2;
    c = P.c_peat + P.c_water*P.f_c * (v+v_new)/2;
    rho = P.density_peat + P.f_rho * P.density_water * P.VMC_0 * (v+v_new)/2;
    % Crank Nicholson approach taking average of times n and n+1
    c_rho = c .* rho;

    %% Update A matrices
    % A1 is dimension 1 update matrix
    % Diagonal elements
    A1_setup(ind_x_inner,3) = -P.alpha_bath/2 - (k(ind_x_inner-1)/2 + k(ind_x_inner) + k(ind_x_inner+1)/2)*fac(1);
    A1_setup(ind_x_upper,3) = -P.alpha_bath/2 - (k(ind_x_upper-1)/2 + k(ind_x_upper)/2 + P.k_sides)*fac(1);
    A1_setup(ind_x_lower,3) = -P.alpha_bath/2 - (k(ind_x_lower+1)/2 + k(ind_x_lower)/2 + P.k_sides)*fac(1);
    % Off-diagonal upper elements
    A1_setup(prod(P.n)+ind_x_inner,3) = (k(ind_x_inner) + k(ind_x_inner+1))/2*fac(1);
    A1_setup(prod(P.n)+ind_x_lower,3) = (k(ind_x_lower) + k(ind_x_lower+1))/2*fac(1);
    % Off-diagonal lower elements
    A1_setup(2*prod(P.n)+ind_x_inner,3) = (k(ind_x_inner) + k(ind_x_inner-1))/2*fac(1);
    A1_setup(2*prod(P.n)+ind_x_upper,3) = (k(ind_x_upper) + k(ind_x_upper-1))/2*fac(1);
    
    % A2 is dimension2 update matrix (to keep this matrix tridiagonal (faster) is operates on the transpose of the T matrix)
    % Diagonal elements
    A2_setup(invtrans(ind_y_inner),3) = -P.alpha_bath/2 - (k(ind_y_inner-P.n(1))/2 + k(ind_y_inner) + k(ind_y_inner+P.n(1))/2)*fac(2);
    A2_setup(invtrans(ind_y_upper),3) = -P.alpha_bath/2 - (k(ind_y_upper-P.n(1))/2 + k(ind_y_upper)/2 + P.k_sides)*fac(2);
    A2_setup(invtrans(ind_y_lower),3) = -P.alpha_bath/2 - (k(ind_y_lower+P.n(1))/2 + k(ind_y_lower)/2 + P.k_sides)*fac(2);
    % For A2_setup we can replace ind_y_upper by trans(ind_x_upper) if region is square
    % Off-diagonal upper elements
    A2_setup(prod(P.n)+invtrans(ind_y_inner),3) = (k(ind_y_inner) + k(ind_y_inner+P.n(1)))/2*fac(2);
    A2_setup(prod(P.n)+invtrans(ind_y_lower),3) = (k(ind_y_lower) + k(ind_y_lower+P.n(1)))/2*fac(2);
    % Off-diagonal lower elements
    A2_setup(2*prod(P.n)+invtrans(ind_y_inner),3) = (k(ind_y_inner) + k(ind_y_inner-P.n(1)))/2*fac(2);
    A2_setup(2*prod(P.n)+invtrans(ind_y_upper),3) = (k(ind_y_upper) + k(ind_y_upper-P.n(1)))/2*fac(2);
    
    % Divide by 2 because of split operator approach
    A1_setup(id1,3) = A1_setup(id1,3) ./ c_rho(A1_setup(id1,1))        * delta_t;
    A2_setup(id2,3) = A2_setup(id2,3) ./ c_rho(trans(A2_setup(id2,1))) * delta_t;
    
    % Create explicit update matrices along dimension 1 and 2
    %    A1_ex = sparse(A1_setup(id1,1), A1_setup(id1,2), A1_setup(id1,3), prod(P.n), prod(P.n));
    %    A2_ex = sparse(A2_setup(id2,1), A2_setup(id2,2), A2_setup(id2,3), prod(P.n), prod(P.n));
    
    if strcmp(computer,'MACI64')
        A1_ex = sponemat + sparse(A1_setup(id1,1),A1_setup(id1,2),A1_setup(id1,3));
        A2_ex = sponemat + sparse(A2_setup(id2,1),A2_setup(id2,2),A2_setup(id2,3));
    else
        A1_ex = sponemat + A1_ex_mask(A1_setup(id1,3));
        A2_ex = sponemat + A2_ex_mask(A2_setup(id2,3));
    end
    
    % Create implicit update matrices along dimension 1 and 2
    A1_imp = 2*sponemat - A1_ex;
    A2_imp = 2*sponemat - A2_ex;
    
    %% Calculate heat sources and sinks
    btmp = (P.alpha_bath*P.Tb + P.h_water*P.VMC_0*P.density_water*dvdt);
    btmp(ind_x_upper) = btmp(ind_x_upper) + P.Ts*P.k_sides*fac(1);
    btmp(ind_x_lower) = btmp(ind_x_lower) + P.Ts*P.k_sides*fac(1);
    btmp(ind_y_upper) = btmp(ind_y_upper) + P.Ts*P.k_sides*fac(2);
    btmp(ind_y_lower) = btmp(ind_y_lower) + P.Ts*P.k_sides*fac(2);
    if nnz(pyro)>0
        btmp(pyro) = btmp(pyro) +  heatsource_pyro(pyro); % [K/cm^3/half time step]
    end
    if nnz(oxy)>0
        btmp(oxy) = btmp(oxy) + heatsource_oxy(oxy); % [K/cm^3/half time step]
    end
    
    if strcmp(computer,'MACI64')
        b1 = sparse(btmp./ c_rho * delta_t);
    else
        b1 = b1_mask(btmp./ c_rho * delta_t);
    end
        
    %% Heat transport
    % First do Dim 1 explicit, dim 2 implicit,
    % then dim2 implicit and dim 1 explicit
    % Transpose dim2 so that A2 matrices can be tri-diagonal (much quicker)
    
    T = sparse(T);    
    
    Ttmp = A1_ex*T + b1;    % Dim 1 explicit
    T = A2_imp \ Ttmp(trans);   % Dim 2 implicit (T transposed)
    
    % End of step 1
    Ttmp = A2_ex * T;           % Dim 2 explicit (T transposed)
    %Ttmp = T(invtrans) + Ttmp(invtrans) + b1;   % Transpose back (can use b1 here in place of b2)
    T = A1_imp \ (Ttmp(invtrans) + b1);          % Dim 1 implicit
    % End of step 2
    
    if t<P.tIgnition/P.dt
        T(ind_igniter) = P.Tigniter;
    end
    
    if P.displayOn % Update the pretty picture of the simulation
        set(t1,'String',['t=' num2str(round(t*P.dt*10)/10,'%05.1f' ) 'mins'])
        set(t2,'String',['Moisture, t=' num2str(round(t*P.dt*10)/10,'%05.1f' ) 'mins'])
        set(t3,'String',['Temperature Distribution, t=' num2str(round(t*P.dt*10)/10,'%05.1f' ) 'mins'])
        set(hp1,'CData', reshape(x,P.n(1),P.n(2)))
        set(hp2,'CData',reshape(P.VMC_0*v_new,P.n(1),P.n(2)))
        set(hp3,'CData',reshape(T+P.abs_zero,P.n(1),P.n(2)))
        drawnow
        %        pause(0.01)
    end
    
    t = t+1;
    v = v_new;  % Update moisture
    
    if (t>P.tIgnition/P.dt && all(T<(50-P.abs_zero)) ) 
        keepGoing=false;
    end  % Stop simulation if no burning cells
end

keep = (1:prod(P.n));

% Save the important times from the simulation
times.pyro = tPyro(keep);
times.char = tChar(keep);
times.burn = tBurn(keep);
times.extinct = tExtinct(keep);

% Trim saved temp and moisture arrays
rec.T = rec.T(1:(t-1),:);
rec.v = rec.v(1:(t-1),:);


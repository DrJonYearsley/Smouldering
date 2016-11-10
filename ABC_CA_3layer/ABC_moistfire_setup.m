function P = ABC_moistfire_setup
% function P = ABC_moistfire_setup
%
% A function that sets all the parameters for a simulation
% All parameters are placed in the structure P
%
% Written Jon Yearsley (jon.yearsley@ucd.ie) July 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P.n=[51, 51]; % width, height Ignition starts along first dimension
P.dt = 0.5; % Time step in mins
P.tMax = 2000;  % Maximum number of iterations [mins]
P.dx = [0.5 0.5]; % Physical step size

P.seed = 0; % if =0 then start with random seed, otherwise use value

P.abs_zero = -273.15; % Absolute zero in degrees C

% Initial conditions ----------
P.Tb = 15 - P.abs_zero; % [K] Temperature of surrounding heat bath
P.Ts = 15 - P.abs_zero; % [K] Temperature of box sides (insulated)

P.Tinit = 15 - P.abs_zero; % Initial temperature of the peat [K]
P.Tigniter = 1000 - P.abs_zero; % Temperature of igniter [K]
P.tIgnition = 30;     % Time igniter is left on [mins]
P.IgniterConfig = 1; % Configuration of the igniter
% =1 Igniter along top side
% =2 Igniter along the middle (dim=1)
% =3 Igniter along the middle (dim=2)
% =4 Igniter a point source in the centre
% =5 Igniter along left side
% =6 Igniter along right side
% =7 Igniter along bottom side

P.VMC_0 = 0.15; % Initial volumetric moisture content [cm^3/cm^3]
% VMC_0=0.04 gives expt velocity of 3.7cm/h = 0.06 cm/min
% VMC_0=0.15 gives expt velocity of 2.5cm/h = 0.04 cm/min

% Parameters for drying
P.drying_model='filkov'; % Can be 'chen' or 'filkov' (default 'filkov')
P.Tmin_evap = 30 - P.abs_zero; % Minimum temp for evaporation
P.VMC_min = 0.0001; % Minimum VMC to consider
P.E_dry_filkov = 42346; %  Activation energy [J/mol] (range=36000-42346)
P.A_dry_filkov = 48200; % Pre-factor [K^0.5 s^-1] (range =48200-450000)
P.E_dry_chen = 3.851E4; %  Activation energy [J/mol]
P.A_dry_chen = exp(8.3); % Pre-factor [s^-1]
P.R = 8.3144621; % Gas constant [J/K/mol]
P.phi = 0.99; % Relative humidity (used in Filkov drying model)

% Smouldering parameters
% Peat char calorific value 5019-6324 kcal/kg max value 8199 kcal/kg (Roy et al 1983)
% (1 cal=4.184 J  therefore ~20 kJ/g, max 2000J/g) Peat calorific
% value (per g) is less (~5000kcal/kg)
% FAO give calorific values for peat of ~20 kJ/g
% (http://www.fao.org/docrep/x5872e/x5872e0b.htm)

% Heat of wood pyrolysis is in range 100-400J/g  (Rath et al 2003)

% Calorific value of peat is ~20kJ/g but not all this may be consummed in
% combustion
P.Hi = 16000; % Effective heat output of oxidising char [J/g peat]
P.Hp = -400; % Heat output of pyrolysis [J/g peat]

% Model parameters for transition probabilities [min and max values]
% Transitions are between states in the burn model:
%         states=fuel, pyrolysis, char, oxidation, ash
P.pyro_Thalf = 150-P.abs_zero; % Temp [K] at which prob of pyrolysis = 50%
P.pyro_A = 10; % Increase in T above Thalf to increase prob to 75%
% If SC_A=0 then pyrlysis starts at precise temperature
P.oxy_Thalf = 400-P.abs_zero;% Temp [K] at which prob of oxidation = 50%
P.oxy_A = 50; % Increase in T above Thalf to increase prob to 75%
P.sigmaI = 0.0001; % Rate of oxidation g peat per min
P.sigmaP = 0.002; % Rate of pyrolysis g peat per min
P.CV = 0.1; % If distribution allows variance to be set...
% Define distribution of pyrolysis and oxidation times
% Can be 'const','exp','unif','gamma','norm'  (const means no stochasticity)
P.pyroDist = 'const';
P.oxyDist = 'const';

% Model Constants ----------

% Thermal conductivity of peat varies from 0.05*10^-2 to 0.15*10^-2 [W/s/K]
P.k_peat = 0.0015; % Thermal conductivity of peat [W/cm/K] = [J/s/cm/K]
P.k_water = 0.006; % Thermal conductivity of peat [W/cm/K]= [J/s/cm/K]
P.k_sides = 0.000; % Thermal conductivity of the sides [W/cm/K]= [J/s/cm/K]
% P.k_air = 0.0024; % (not used) Thermal conductivity of air [W/cm/K]= [J/s/cm/K]
% 0.001 W/cm/K is conductivity of brick/wood/gypsum board

P.c_peat = 0.92;   % Specific heat of peat [J/g/K] (Oke 1997, Benscoter et al 2011)
P.c_water = 4.187; % Specific heat of water [J/g/K]
P.c_sides = 1.0;   % Specific heat of the sides [J/g/K]
% Specific heat of plaster board ~0.72 J/cm^3/K, density~0.8g/cm^3

P.h_water = 2270; % Latent heat of evaporation [J/g]
P.h_ice = 334; % Latent heat of fusion [J/g]

% Moisture dependent heat capacities and conductances
% Conductivity is k = k_peat + f_k * k_water * VMC
% Heat Capacity is c = c_peat + f_c * c_water * VMC
% Desnity is rho = rho_peat + f_rho * rho_water * VMC
% Dissanayaka et al (2011) Vadose Zone Journal doi:10.2136/vzj2011.0092
P.f_k = 0.89;
P.f_c = 0.88;
P.f_rho = 1;

% P.f_k = 0;
% P.f_c = 0;
% P.f_rho = 0;
% P.h_water = 0;

% Densities
% Peat = 151 +/- 3 kg /m^3 [10^-3 g/cm^3]
% Ash  = 36 +/- 2 kg /m^3  [10^-3 g/cm^3]
% Char = 189 +/- 4 kg /m^3 [10^-3 g/cm^3]
P.density_peat = 151 *10^-3; % g cm^-3
P.density_water = 1;         % g cm^-3   [or 1000 kg/m^3]

% A factor governing heat losses throughout the sample to the heat bath
P.alpha_bath = 0.08 * P.density_peat*P.c_peat/60;  % Units [J/ cm^3 / s]
% Cooling rate from Nuria's burns is approximately 0.003 K/min for late
% cooling and 0.04 K/min just after peak

% References:
% Roy et al 1983 The pyrloysis of peat, Journal of analytical and applied
%    pyrolysis, Vol 5, p 261-332
% Dissanayaka et al (2011) Vadose Zone Journal doi:10.2136/vzj2011.0092
% Rath et al 2003 Heat of wood pyrloysis, Fuel, Vol 82, p81-91
% Benscoter et al (2011), Interactive effects of vegetation, soil moisture
%    and bulk density on depth of burning of thick organic soils
%    Int J of Wildland Fire Vol 20, p 418-429

% Select points to record temperature profile
P.tRec = [10,50,100,200];
%P.indRec = floor(P.n(2)/2)*P.n(1)+floor(linspace(1,P.n(1),10));
P.indRec = floor(P.n(2)/2)*P.n(1)+floor(linspace(1,P.n(1),10));
[i,j]=ind2sub(P.n,P.indRec);
P.xRec = (diag(P.dx)*[i;j])';

%% Set random number generator seed
if P.seed==0
    rng('shuffle')
    P.seed = rng;
else
    rng(P.seed)
    P.seed = rng;
end
function epsilon = ABC_distance(result, target)
% FUNCTION espilon = ABC_distance(result, target)
% 
% A function to calculate a distance measure between the
% experimental data and the results from one simulation
%
% Inputs: 
%    result   The simulation data
%    target   A structure containing the experimental data
%      target.tol   = the acceptable tolerance for each observation
%      target.point = the data to be compared to the simulation data
% 
% Output:
%    epsilon The distance measure. In this case it is the Euclidean
%            distance after rescaling by target.tol 
%
% Written: Jon Yearsley July 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon  = sqrt(sum((result(:)-target.point(:)).^2./target.tol(:).^2));

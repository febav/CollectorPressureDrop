function rho=densityGlyMixAndWat(x,temp)

% function rho=densityGlyMixAndWat(x,temp)
% This function computes the density of water and glycol/water mixtures
% given glycol content [%] and fluid temperature [degC] as input.
% Note that the temperature can be also passed as vector. In this case
% the function output will be a vector.

if x>0   % density for glycol/water mixtures, DTU KEM (35%<x<50%)
    rho=1013-0.2682*temp+0.7225*x-0.00194*temp.^2-0.004964*temp*x; % [kg/m3]
elseif x==0     % if glycol%=0, switch to water equation (Simon's equation)
    rho=1000.6-0.0128*temp.^1.76; % [kg/m3]
else
    error('Glycol % must be >=0!');
end
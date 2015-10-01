function nu=viscosityGlyMixAndWat(x,temp)

% function nu=viscosityGlyMixAndWat(x,temp)
% This function computes the viscosity of water and glycol/water mixtures
% given glycol content [%] and fluid temperature [degC] as input.
% Note that the temperature can be also passed as vector. In this case
% the function output will be a vector.

if x>35   % viscosity for glycol/water mixtures, DTU KEM (=exp mu/exp rho)
    nu=(-2.881-6.721*10^-3*temp+0.2839*x+1.959*10^-3*temp.^2+...
        -7.036*10^-3*x*temp-1.883*10^-5*temp.^3+...
        4.862*10^-5*x*temp.^2)/1000./...
        densityGlyMixAndWat(x,temp); % [m2/s]
elseif x==0  % if glycol%=0, switch to water equation (=Kestin nu/Simon's rho)
    nu=1.002/1000*10.^((20-temp)./(temp+96).*(1.2378-1.303*10^-3*(20-temp)+...
        3.06*10^-6*(20-temp).^2+2.55*10^-8*(20-temp).^3))./...
        densityGlyMixAndWat(x,temp); % [m2/s]
else
    error('Glycol % is not valid!');
end
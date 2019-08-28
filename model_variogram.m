function V=model_variogram(h,c,a,n,type)
% variogram using the bounded linear and spherical models
% HPM 10/27/09
% INPUT: h = lags at which modeled estimates are made
%        c = variogram sill = sig^2 (model parameter)
%        a = variogram range (model parameter)
%        n = nugget (model parameter)
%      type = 'L' for linear, 'S' for spherical

% first calculate variogram for lags less than range
ind=find(h<=a);
if strcmp(type,'L')
    V(ind)=c*h(ind)/a; % bounded linear
elseif strcmp(type,'S')
    V(ind)=c*(3*h(ind)/(2*a)-1/2*(h(ind)/a).^3); % spherical
elseif strcmp(type,'E')
    V=c*( 1-exp(-h/a));
elseif strcmp(type,'P')
    V(ind)=w*(1-cos((2*pi*h)/omega));
end    


% now define variogram for lags greater than range
if ~strcmp(type,'E')
 ind2=find(h>a); % find points greater than range
V(ind2)=c; % set equal to sill
end
V=V+n; % add nugget

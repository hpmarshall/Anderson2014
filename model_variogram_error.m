function chi2w=model_variogram_error(h,V,w,c,a,n,type)
% variogram using the bounded linear and spherical models
% HPM 10/27/09
% INPUT: h = lags at which modeled estimates are made
%            V = experimental variogram
%       w = weights for each point
%        c = variogram sill = sig^2 (model parameter)
%        a = variogram range (model parameter)
%        n = nugget (model parameter)
%      type = 'L' for linear, 'S' for spherical

% first calculate variogram for lags less than range
ind=find(h<=a);
if strcmp(type,'L')
    Vmod(ind)=c*h(ind)/a; % bounded linear

else
    Vmod(ind)=c*(3*h(ind)/(2*a)-1/2*(h(ind)/a).^3); % spherical
end

% now define variogram for lags greater than range
ind2=find(h>a); % find points greater than range
Vmod(ind2)=c; % set equal to sill
Vmod=Vmod+n; % add nugget
chi2w=sqrt(mean(w/sum(w).*(V(:)-Vmod(:)).^2));
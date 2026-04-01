function f_daugment(path,spec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes as given the nonlinear parameters of an RCNL model,
% obtains the mean consumer valuation (delta) using the contraction 
% mapping, and then saves it along other statistics in a new data file.
%
% Called by:
%   - main_v6.m
% Calls the following user-specified functions:
%   - main_data.m 
%   - rcnl_meanval.m
% Loads the following data:
%   - dres_gmm2 (file from MW [2017])
% Saved data used in:
%   - main_supply_bind.m
%   - main_supply_nonbind.m
%   - impute_bertrand.m
%   - main_supply_bind_nopool.m
%   - main_supply_bind_x.m
%   - f_getdev.m
%   - f_pi_m
%   - f_impute_mc.m
%   - results_baseanalysis.m
%   - results_cfmergers.m
%   - results_cfanalysis.m
%   - idgraph.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Loading demand results
cd(strcat(path.data1))
%load dres_gmm2
load small_dresgmm2
derMat=derMat_2;
elasMat=elasMat_2;
theta2=theta2_2;
alpha=theta1_2(1);
rho=rho_2;
if strcmp(spec.dfolder,'RCNL1') || strcmp(spec.dfolder,'RCNL2')
    theti = [1;2;3];
    thetj = [2;2;2];
elseif strcmp(spec.dfolder,'RCNL3') || strcmp(spec.dfolder,'RCNL4')
    theti = [2;3;4;5];
    thetj = [2;2;2;2];
end
theta2w = full(sparse(theti,thetj,theta2));
clear k temp
cd(strcat(path.code1))

% Structure for demand results
daugfile.rho = rho;
daugfile.alpha = alpha;
daugfile.derMat = derMat;
daugfile.elasMat = elasMat;
if strcmp(spec.dfolder,'NL1')==1
    daugfile.theta2w = 0;
    daugfile.theti = 0;
    daugfile.thetj = 0;
    daugfile.theta2 = 0;
else
    daugfile.theta2w = theta2w;
    daugfile.theti = theti;
    daugfile.thetj = thetj;
    daugfile.theta2 = theta2;
end

% Grabbing data
[vars,idmatrix] = main_data(path,spec);

% Contraction mapping
[delta,mu,ai] = rcnl_meanval(vars,idmatrix,daugfile);
pcoefi = alpha+ai;

% Mean non-price valuation
deltanp = delta - alpha*vars.p_jt;

% Fixed effects coefficients
dfecoef = (vars.fesd'*vars.fesd)\vars.fesd'*deltanp;
dprodfecoef = dfecoef(1:length(unique(idmatrix.prodid)));
ddatefecoef = [0;dfecoef(1+length(unique(idmatrix.prodid)):end)];

% Residual (mean unobserved quality) and city-specific component
xi = deltanp - vars.fesd*dfecoef;
citefe = cr_dum(idmatrix.cityid);
dcityfecoef = (citefe'*citefe)\citefe'*xi;

% Structure for demand results
daugfile.delta = delta;
daugfile.mu = mu;
daugfile.ai = ai;
daugfile.pcoefi = pcoefi;
daugfile.deltanp = deltanp;
daugfile.dprodfecoef = dprodfecoef;
daugfile.ddatefecoef = ddatefecoef;
daugfile.dcityfecoef = dcityfecoef;
daugfile.xi = xi;

% Saving augmented demand results
cd(strcat(path.data2))
save daugfile daugfile
cd(path.code1)



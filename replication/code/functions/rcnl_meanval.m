function [f,g,h] = rcnl_meanval(vars,idmatrix,daugment)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is analogous to the Nevo "meanval" function.  The
% contraction mapping is outsourced to a C program, which 
% increases speed in part by looping efficiently over markets.  This
% function returns the vector of mean valuations (delta) and also a matrix
% of consumer-specific deviations (mu).  The latter object is used in
% supply-side moments and counter-factuals.
%
% Called by:
%   - f_daugment.m
%
% Calls the following user-specified functions:
%   - f_mu.m 
%   - contrMap_rcnl.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Starting value for contraction mapping
temp = log(vars.s_jt)-log(vars.outshr)-daugment.rho*vars.logcondshr;
expmvalold = exp(temp);

[mu,ai] = f_mu(daugment.theta2w,vars.ns,vars.dfull,vars.x2);
expmu = exp(mu);

% Contraction mapping is done in C 
maxi = 2500;
tol  = 1e-14;
[expmval,b,~] = contrMap_rcnl(tol,maxi,expmvalold,expmu,...
    idmatrix.cdindex,vars.s_jt,daugment.rho);

% Warning message if contraction mapping fails
if max(b)>tol || sum(isnan(expmval))>0
    disp('WARNING: CONTRACTION MAPPING CONVERGENCE CRITERION FAILS');
end
    
f = log(expmval);
g = mu;
h = ai;



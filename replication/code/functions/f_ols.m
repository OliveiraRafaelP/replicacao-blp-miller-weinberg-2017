function [gamma,se,var,tstat] = f_ols(y,X,cv)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function does OLS with a y vector, X matrix, and cv to define
% clusters for the SEs.
%
% Called by:
%   - results_costregs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
% OLS regression with clustered standard errors
gamma = inv(X'*X)*X'*y;
resid = y - X*gamma;
[~,ncol] = size(X);
sqerr = zeros(ncol,ncol);
for m=1:max(cv)
    sqerr=sqerr+X(cv==m,:)'*(resid(cv==m)*resid(cv==m)')*X(cv==m,:);
end
var = inv(X'*X)*sqerr*inv(X'*X);
se  = sqrt(diag(var));
tstat = gamma ./ se;

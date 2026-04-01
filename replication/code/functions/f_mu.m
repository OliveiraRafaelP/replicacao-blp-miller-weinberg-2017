function [f, g] = f_mu(theta2w,ns,dfull,x2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is a renamed version of the Nevo function "mufunc".  It
% has been edited to return the individual-specific component of
% indirect utility (mu) in matrix form and the individual-specific
% component of the price coefficient (ai) also in matrix form.  
%
% Called by:
%   - f_meanval.m
%   - cf_foc_partial.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,k] = size(x2);
j = size(theta2w,2)-1;
mu = zeros(n,ns);
ai = zeros(n,ns);

for i = 1:ns
    d_i = dfull(:,i:ns:j*ns);
    mu(:,i) = x2.*(d_i*theta2w(:,2:j+1)')*ones(k,1);
    ai(:,i) = d_i*theta2w(1,2:j+1)';
end

f = mu;
g = ai;



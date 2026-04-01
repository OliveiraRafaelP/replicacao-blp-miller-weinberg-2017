function [zero,svec,der1,bMarkup] = cf_foc_partial(p,alpha,rho,deltanp2,cost,owner,pfix,...
    popt,x2M,pcoefiM,dfullM,theta2w,ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the value of the first order conditions in a
% market, given a set of conduct/demand/cost conditions.  The primary use
% is in the computation of equilibrium.  Also used to obtain fringe
% markups.
%
% Called by at least:
%   - f_impute_mc.m
%   - f_pi_m.m
%   - f_impute_mc.m
% Calls the following user-specified functions:
%   - f_mu.m
%   - rcnl_indsh.m
%   - rcnl_der1.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Market prices
pvec = zeros(size(popt));
pvec(popt==1) = p;
pvec(popt==0) = pfix;

% Adjusted mean valuations 
delta2 = deltanp2 + alpha*pvec;

% Next steps demand on whether NL (faster) or RCNL (slower); NL has ns==0
if ns==0
    
    % Market shares
    Dg = sum  ( exp(delta2 / (1-rho)));
    scond = exp(delta2/(1-rho)) / Dg;
    sgrup = Dg^(1-rho) / (1 + Dg^(1-rho));
    svec = scond*sgrup;
    
    % Derivatives
    cross = -alpha*(svec + rho/(1-rho)*scond)*svec';
    own = alpha /(1-rho) * svec .*( 1 - rho*scond - (1-rho)*svec);
    der1 = cross - diag(diag(cross)) + diag(own);

else
    
    %Adjusting consumer-specific preference deviations
    x2M(:,1) = pvec;
    [mu2,ai] = f_mu(theta2w,ns,dfullM,x2M);
    
    
    % Consumer-specific choice probabilities and market shares
    [sharei,scondi,sgroupi]=rcnl_indsh(exp(delta2),exp(mu2),rho,0,0,1);
    svec = mean(sharei,2);
    
    % Derivatives
    der1 = rcnl_der1(pcoefiM,sharei,scondi,sgroupi,rho);
    
end

%Obtain implied brewer markup
bMarkup   = - (owner.*der1(popt==1,popt==1)')\svec(popt==1);

%First order condition
zero = p - cost - bMarkup;





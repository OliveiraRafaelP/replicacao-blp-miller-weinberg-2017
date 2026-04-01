function [c, ceq]= icc(inc, p_n, coal, firmid, xi, mc, alpha, delta, m, outside, prodid, optsol)

% Description: calculates value of the ICC at the firm level
% Inputs:
    % Let J be the number of products (excluding the outside good)
    % inc: scalar increase in price above Nash
    % p_n: (Jx1) vector of Nash equilibrium prices
    % coal: (number of firms x 1) dummy vector, equals 1 if that firm is in 
    % the coalition
    % firmid: (number of firms x 1) vector of integer firm ID numbers 
    % xi: (Jx1) vector of quality parameters
    % mc: (Jx1) vector of marginal costs
    % alpha: scalar price coefficient
    % delta: scalar discount factor
    % m: scalar market size
    % outside: scalar normalization for outside good value
    % prodid: (Jx1) vector that identifies which firm goes with which
    % product
    % optsol: options for the fsolve routine
    
% Output:
    % c: (number of coalition firms x 1) vector of ongoing value of
    %  cheating and punishment less the ongoing value of collusion
    % ceq: empty place holder matrix

% Profits if stay in the coalition
[~, pi_c] = pipinc(inc, p_n, coal, xi, mc, alpha, m, outside, prodid, optsol);
pi_c = pi_c(coal == 1);

% Nash profits
[~, ~, pi_n]= logit(p_n, xi, mc, alpha, m, outside, prodid);
pi_n = pi_n(coal == 1);

% Deviation profits (firm-by-firm)
coalfirms = firmid(coal == 1);
pi_d = zeros(size(coalfirms));
for i = 1:length(coalfirms)
    coaltemp = coal;
    coaltemp(coalfirms(i)) = 0;
    [~, pi_temp] = pipinc(inc, p_n, coaltemp, xi, mc, alpha, m, outside, prodid, optsol);
    pi_d(i) = pi_temp(coalfirms(i));
end

% Form ICC for inequality constraints
% To satisfy the ICCs, these should be less than or equal to zero
c = (pi_d + (delta/(1-delta)) * pi_n) - (1/(1-delta)) * pi_c;

% Equality constraints placeholder
ceq = [];




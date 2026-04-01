function [loss, p_n, mc]= recovery_1(inchat, p_c, s_c, pi_c, coal, firmid, xi, alpha, delta, m, outside, prodid, optsol)

% Description: used to recover the supermarkup with binding ICC
% Inputs:
    % Let J be the number of products (excluding the outside good)
    % inchat: candidate scalar increase in price above Nash
    % p_c: (Jx1) vector of PLE equilibrium prices
    % coal: (number of firms x 1) dummy vector, equals 1 if that firm is in 
    % the coalition
    % firmid: (number of firms x 1) vector of integer firm ID numbers 
    % xi: (Jx1) vector of quality parameters
    % alpha: scalar price coefficient
    % delta: scalar discount factor
    % m: scalar market size
    % outside: scalar normalization for outside good value
    % prodid: (Jx1) vector that identifies which firm goes with which
    % product
    % optsol: options for the fsolve routine
    
% Output:
    % loss function
    % marginal costs
    % prices in nash

    
% Profits if stay in the coalition
pi_c = pi_c(coal == 1);

% Nash prices and marginal costs ---- USES THE MC IDENTIFICATION PROCEDURE
pfix = p_c - inchat;
p_n = zeros(size(p_c));
mc = zeros(size(p_c));
if min(coal) == 0
    owner = f_ownMat(firmid(coal==0));
    dqdp = s_c(coal==0)*s_c(coal==0)';
    dqdp = dqdp - diag(diag(dqdp)) - diag(s_c(coal==0).*(1-s_c(coal==0)));
    dqdp = -dqdp*alpha;
    mc(coal==0) = p_c(coal==0) + inv(owner.*dqdp')*s_c(coal==0);
    clear owner dqdp
    f = @(p)nashfoc(p, pfix, coal, xi, mc, alpha, m, outside, prodid);
   [p_fringe] = fsolve(f, mc(coal == 0), optsol);
    p_n(coal == 1) = pfix(coal == 1);
    p_n(coal == 0) = p_fringe;
    [s_n, ~, ~]= logit(p_n, xi, mc, alpha, m, outside, prodid);
    owner = f_ownMat(firmid(coal==1));
    dqdp = s_n(coal==1)*s_n(coal==1)';
    dqdp = dqdp - diag(diag(dqdp)) - diag(s_n(coal==1).*(1-s_n(coal==1)));
    dqdp = -dqdp*alpha;
    mc(coal==1) = p_n(coal==1) + inv(owner.*dqdp')*s_n(coal==1);
else
    p_n = pfix;
    [s_n, ~, ~]= logit(p_n, xi, mc, alpha, m, outside, prodid);    
    owner = f_ownMat(firmid(coal==1));
    dqdp = s_n(coal==1)*s_n(coal==1)';
    dqdp = dqdp - diag(diag(dqdp)) - diag(s_n(coal==1).*(1-s_n(coal==1)));
    dqdp = -dqdp*alpha;
    mc(coal==1) = p_n(coal==1) + inv(owner.*dqdp')*s_n(coal==1);    
end
[~, ~, pi_n]= logit(p_n, xi, mc, alpha, m, outside, prodid);
pi_n = pi_n(coal == 1);


% Deviation profits (firm-by-firm)
coalfirms = firmid(coal == 1);
pi_d = zeros(size(coalfirms));
for i = 1:length(coalfirms)
    coaltemp = coal;
    coaltemp(coalfirms(i)) = 0;
    [~, pi_temp] = pipinc(inchat, p_n, coaltemp, xi, mc, alpha, m, outside, prodid, optsol);
    pi_d(i) = pi_temp(coalfirms(i));
end

% Form ICC for inequality constraints
% To satisfy the ICCs, these should be less than or equal to zero
g = (pi_d + (delta/(1-delta)) * pi_n) - (1/(1-delta)) * pi_c;

% Loss equals 0 if ICC binds
loss = (max(g)^2);


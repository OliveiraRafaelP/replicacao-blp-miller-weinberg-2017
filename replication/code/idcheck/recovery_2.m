function [loss, p_n, mc]= recovery_2(inchat, p_c, s_c, coal, firmid, leadid, xi, alpha, m, outside, prodid, optsol)

% Description: used to recover the supermarkup with nonbinding ICC
% Inputs:
    % Let J be the number of products (excluding the outside good)
    % inchat: candidate scalar increase in price above Nash
    % p_c: (Jx1) vector of PLE equilibrium prices
    % coal: (number of firms x 1) dummy vector, equals 1 if that firm is in 
    % the coalition
    % firmid: (number of firms x 1) vector of integer firm ID numbers 
    % xi: (Jx1) vector of quality parameters
    % alpha: scalar price coefficient
    % m: scalar market size
    % outside: scalar normalization for outside good value
    % prodid: (Jx1) vector that identifies which firm goes with which
    % product
    % optsol: options for the fsolve routine
    
% Output:
    % loss function
    % marginal costs
    % prices in nash


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


%Getting dpi_c / dm via numerical derivative

eps = 1e-4;

inchatx = inchat - eps/2;
pfixx = p_n + inchatx;
px_c = zeros(size(pfixx));
f = @(p)nashfoc(p, pfixx, coal, xi, mc, alpha, m, outside, prodid);
[px_fringe] = fsolve(f, mc(coal == 0), optsol);
px_c(coal == 1) = pfixx(coal == 1);
px_c(coal == 0) = px_fringe;
[~, ~, pix_c1]= logit(px_c, xi, mc, alpha, m, outside, prodid);
clear inchatx pfixx px_c px_fringe px_c

inchatx = inchat + eps;
pfixx = p_n + inchatx;
px_c = zeros(size(pfixx));
f = @(p)nashfoc(p, pfixx, coal, xi, mc, alpha, m, outside, prodid);
[px_fringe] = fsolve(f, mc(coal == 0), optsol);
px_c(coal == 1) = pfixx(coal == 1);
px_c(coal == 0) = px_fringe;
[~, ~, pix_c2]= logit(px_c, xi, mc, alpha, m, outside, prodid);
clear inchatx pfixx px_c px_fringe px_c

dpidm = sum(pix_c2(leadid==1)-pix_c1(leadid==1)) / eps;

% Loss equals 0 if profit is flat in supermarkup
loss = dpidm^2;


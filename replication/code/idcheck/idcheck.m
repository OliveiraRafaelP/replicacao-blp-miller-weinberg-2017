function [mc1,mc2,smarkup1,smarkup2,bind] = idcheck(sc,sout,mark,coalid,delta,alt)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This code tests identification given observables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Obtain structural parameters using a logit/Bertrand model. We
% start with shares and a markup (prices normalized to one) in order to
% keep the structural parameters in a plausible range. This calibration
% methodology tracks MRRS (2016,2017).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prices
p0 = ones(size(sc));

%base calibration -- results structure
calbase = cal_logit(sc,sout,mark,p0);

% alternative calibration using draws on quality, cost, not the input data
nprods = length(sc);
calalt.alpha = -1;
calalt.xi = 1 + rand(nprods,1);
calalt.mc = rand(nprods,1);

pfix = zeros(nprods,1);
optsol = optimoptions(@fsolve,'Disp','off');
f = @(p)nashfoc(p, pfix, pfix, calalt.xi, calalt.mc, calalt.alpha, 1, 0, [1:nprods]');
calalt.p_n = fsolve(f, calalt.mc, optsol);

calalt.mark = calalt.p_n - calalt.mc;
num = exp(calalt.xi + calalt.alpha*calalt.p_n);
denom = sum(num) + 1;
calalt.s_n = num / denom;
calalt.pi = calalt.mark .* calalt.s_n;

% selecting the approach to calibration
if alt==1
    cal = calalt;
else
    cal = calbase;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Simulate the PLE given the structural parameters. The resulting
% prices and shares will be treated as observables in the next step. PLE
% simulation requires an assumption on the discount factor, firm
% identities, and the coalition structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Struction on DGP
prodid = [1:length(coalid)]';
leadid = prodid == 1;

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Bounds
lb = -5;
ub = 10 * max(cal.p_n) - max(cal.p_n);

% Starting value
start = 1.5 * max(cal.p_n) - max(cal.p_n);

%%% Optimization options
optsol = optimoptions(@fsolve,'Disp','off');
optcon = optimoptions(@fmincon,'Disp','off');

% Collect constraint function and function to be maximized
g = @(incx)icc(incx, cal.p_n, coalid, prodid, cal.xi, cal.mc, cal.alpha, delta, 1, 0, prodid, optsol);
fnx = @(inc)leadpipinc(inc, cal.p_n, coalid, leadid, cal.xi, cal.mc, cal.alpha, 1, 0, prodid, optsol);

% Do minimization and record output
[inc_c, ~, flag_c, ~, lambda] = fmincon(fnx, start, A, b, Aeq, beq, lb, ub, g, optcon);
[p_c, pi_c] = pipinc(inc_c, cal.p_n, coalid, cal.xi, cal.mc, cal.alpha, 1, 0, prodid, optsol);
[s_c, ~] = logit(p_c, cal.xi, cal.mc, cal.alpha, 1, 0, prodid);
lagrangian = lambda.ineqnonlin;

% Structure with results of DGP
dgp.p_c = p_c;
dgp.s_c = s_c;
dgp.pi_c = pi_c;
dgp.inc_c = inc_c;
dgp.lagrangian = lagrangian;
dgp.coalid = coalid;
dgp.prodid = prodid;
dgp.leadid = leadid;
dgp.delta = delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Recover the supermarkup (inc_c) given observed data
%   - only do this if we converge on a PLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if flag_c <= 0 
    
    disp('CONVERGENCE PROBLEM IN DGP') ; 
    mc1 = [];
    mc2 = [];
    smarkup1=[];
    smarkup2=[];
    bind=[];
    
else 
    
    if max(lagrangian)>1e-3
        
        bind=1;
        fs = @(inchat)recovery_1(inchat, dgp.p_c, dgp.s_c, dgp.pi_c, dgp.coalid, dgp.prodid, cal.xi, cal.alpha, dgp.delta, 1, 0, dgp.prodid, optsol);
        [smarkup,~,~,~] = fminsearch(fs,0.01);
        [~, ~, mc] = recovery_1(smarkup, dgp.p_c, dgp.s_c, dgp.pi_c, dgp.coalid, dgp.prodid, cal.xi, cal.alpha, dgp.delta, 1, 0, dgp.prodid, optsol);
        
    else 
        
        bind=0;
        fs = @(inchat)recovery_2(inchat, dgp.p_c, dgp.s_c, dgp.coalid, dgp.prodid, dgp.leadid, cal.xi, cal.alpha, 1, 0, dgp.prodid, optsol);
        [smarkup,~,~,~] = fminsearch(fs,0.01);
        [~, ~, mc] = recovery_2(smarkup, dgp.p_c, dgp.s_c, dgp.coalid, dgp.prodid, dgp.leadid, cal.xi, cal.alpha, 1, 0, dgp.prodid, optsol);
        
    end
   
    mc1 = cal.mc;
    mc2 = mc;
    smarkup1=dgp.inc_c;
    smarkup2=smarkup;
    
end

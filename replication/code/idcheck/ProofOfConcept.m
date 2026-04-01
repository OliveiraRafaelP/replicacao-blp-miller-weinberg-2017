%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This code conducts an exercise on identification using PLE with 
%%%% logit demand with (i) one region, and then (ii) multiple regions. The
%%%% steps include:
%%%%  (1) obtain structural parameters using a logit/Bertand calibration
%%%%  (2) simulate PLE as a constrained maximization problem
%%%%  (3) Take the PLE data and recover the supermarkup(s) and costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set path
cd(strcat(path.code,'/','idcheck'))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Obtain structural parameters using a logit/Bertrand model. We
% start with shares and a markup (prices normalized to one) in order to
% keep the structural parameters in a plausible range. This calibration
% methodology tracks MRRS (2016,2017).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%conditional shares sum to one
sc0 = [0.4 0.2 0.2 0.1 0.1]';
sc0 = sc0/sum(sc0);

%outside good hsare
sout0 = 0.2;

%markup of first product
m0 = 0.5;

%prices
p0 = ones(size(sc0));

%calibration results structure
cal = cal_logit(sc0,sout0,m0,p0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Simulate the PLE given the structural parameters. The resulting
% prices and shares will be treated as observables in the next step. PLE
% simulation requires an assumption on the discount factor, firm
% identities, and the coalition structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Struction on DGP
delta = 0.8;
coalid = [1 1 1 0 0]';
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
[inc_c, fval_c, flag_c, ~, lambda] = fmincon(fnx, start, A, b, Aeq, beq, lb, ub, g, optcon);
if flag_c <= 0 ; disp('CONVERGENCE PROBLEM IN DGP') ; end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if max(lagrangian)>1e-4
    
    disp('ICC Binds: Using Recovery 1');
    fs = @(inchat)recovery_1(inchat, dgp.p_c, dgp.s_c, dgp.pi_c, dgp.coalid, dgp.prodid, cal.xi, cal.alpha, dgp.delta, 1, 0, dgp.prodid, optsol);
    [smarkup,sfvala,~,~] = fminsearch(fs,0.01);
    [~, p_n, mc] = recovery_1(smarkup, dgp.p_c, dgp.s_c, dgp.pi_c, dgp.coalid, dgp.prodid, cal.xi, cal.alpha, dgp.delta, 1, 0, dgp.prodid, optsol);

else
    
    disp('ICC Does Not Bind: Using Recovery 2');
    fs = @(inchat)recovery_2(inchat, dgp.p_c, dgp.s_c, dgp.coalid, dgp.prodid, dgp.leadid, cal.xi, cal.alpha, 1, 0, dgp.prodid, optsol);
    [smarkup,sfvala,~,~] = fminsearch(fs,0.01);
    [~, p_n, mc] = recovery_2(smarkup, dgp.p_c, dgp.s_c, dgp.coalid, dgp.prodid, dgp.leadid, cal.xi, cal.alpha, 1, 0, dgp.prodid, optsol);
    
end

impute.p_n = p_n;
impute.mc = mc;
impute.smarkup = smarkup;

[cal.mc impute.mc]
[dgp.inc_c impute.smarkup]



cd(path.code)


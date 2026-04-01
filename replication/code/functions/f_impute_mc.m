function [mc,pnash,snash] = ...
    f_impute_mc(sm,mktid,vars,ids,daugfile,spec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function imputes marginal costs for some supermarkup, m, for a
% given period*city "market" defined by cdid==m. It also imputes Bertrand
% prices, shares, and profit, as well as PLE profit.
%
% Called by:
%   - f_loss_nonbind.m
%   - f_loss_bind.m
%   - impute_bertrand.m
%   - f_loss_bind_x.m
% Calls the following user-specified functions:
%   - f_ownMat.m
%   - cf_foc_partial.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm and coalition identity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdid = ids.cdid;
firmid = ids.firmid(cdid==mktid);
sizeid = ids.sizeid(cdid==mktid);
coalid = ids.coalid(cdid==mktid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracting relevant prices, shell files, other data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_jt = vars.p_jt(cdid==mktid);
mc = zeros(size(p_jt));
pnash = zeros(size(p_jt));


pcoefi = daugfile.pcoefi(ids.obsindemand==1,:);
deltanp = daugfile.deltanp(ids.obsindemand==1);

x2M = vars.x2(cdid==mktid,:);
pcoefiM = pcoefi(cdid==mktid,:);
deltanpM = deltanp(cdid==mktid);
dfullM = vars.dfull(cdid==mktid,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first step is to recover the marginal costs of the fringe firms,
% which is done on the basis of the first order conditions.
%   - shorter code is possible but this minimizes duplicates code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fringe markups and marginal costs
owner   = f_ownMat(firmid(coalid==0));
p = p_jt(coalid==0);
popt = 1- coalid;
pfix = p_jt(coalid==1);
[~,~,~,bvec] = cf_foc_partial(p,daugfile.alpha,daugfile.rho,deltanpM,zeros(size(p)),...
    owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns);
mc(coalid==0) = p_jt(coalid==0) - bvec;
clear owner p pot pfix bvec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The second step is to recover Nash prices.  This is trivial for coalition
% firms and pre-coalition markets; computation is required for fringe
% firms in coalition markets.  The computation is the main bottleneck.
%   - appears robust to starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bertrand prices for coalition firms
if spec.bysize==0
    pnash(coalid==1) = p_jt(coalid==1) - sm;
elseif spec.bysize==1
    pnash(coalid==1 & sizeid==1) = p_jt(coalid==1 & sizeid==1) - sm(1);
    pnash(coalid==1 & sizeid==2) = p_jt(coalid==1 & sizeid==2) - sm(2);
    pnash(coalid==1 & sizeid==3) = p_jt(coalid==1 & sizeid==3) - sm(3);
elseif spec.bysize==2
    pnash(coalid==1 & sizeid==1) = p_jt(coalid==1 & sizeid==1) - sm(1);
    pnash(coalid==1 & sizeid==2) = p_jt(coalid==1 & sizeid==2) - sm(1);
    pnash(coalid==1 & sizeid==3) = p_jt(coalid==1 & sizeid==3) - sm(2);
end

% Bertrand prices for fringe firms
owner   = f_ownMat(firmid(coalid==0));
p = p_jt(coalid==0);
popt = 1- coalid;
pfix = pnash(coalid==1);
f = @(p)cf_foc_partial(p,daugfile.alpha,daugfile.rho,deltanpM,mc(coalid==0),...
    owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns);
[fpn,~,test,~] = fsolve(f,p,spec.options.cf);
pnash(coalid==0) = fpn;
clear owner p popt pfix fpn


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The third step is to calculate the shares and demand derivatives 
% evaluated at Nash prices, and back out the implied marginal costs for 
% coalition firms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nash shares and derivatives
owner   = f_ownMat(firmid);
p = pnash;
popt = ones(size(p));
pfix = [];
[~,snash,der1,~] = cf_foc_partial(p,daugfile.alpha,daugfile.rho,deltanpM,mc,...
    owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns);      
clear owner p popt pfix

% Nash markups and marginal costs for coalition firms
der2  = der1(coalid==1,coalid==1);
owner = f_ownMat(firmid(coalid==1));
bvec  = - (owner.*(der2'))\snash(coalid==1);
mc(coalid==1) = pnash(coalid==1) - bvec;
clear der1 der2 owner bvec




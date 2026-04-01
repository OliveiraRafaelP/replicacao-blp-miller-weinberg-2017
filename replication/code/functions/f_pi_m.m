function [price,share,pi,pdev,sdev,pidev] = ...
    f_pi_m(sm,mktid,pnashx,mcx,vars,ids,daugfile,spec,getdev,restr,modm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the vector of prices, shares, and profits realized 
% with a given supermarkup, sm, given Bertrand prices, marginal costs, and 
% demand, and for a specific period*city "market" definded by cdid==mktid.
% If requested via getdev==1 then it also obtains the corresponding
% deviation objects.
%
% Called by:
%   - f_loss_nonbind.m
%   - f_loss_bind.m
%   - f_icc_eval.m
%   - f_loss_bind_x.m
% Calls the following user-specified functions:
%   - f_ownMat.m
%   - cf_foc_partial.m
%   - f_getdev.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm and coalition identity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdid = ids.cdid;
firmid = ids.firmid(cdid==mktid);
sizeid = ids.sizeid(cdid==mktid);
coalid = ids.coalid(cdid==mktid);
brndid = ids.brndid(cdid==mktid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shell files, demand data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

price = zeros(size(firmid));

pcoefi = daugfile.pcoefi(ids.obsindemand==1,:);
deltanp = daugfile.deltanp(ids.obsindemand==1);

x2M = vars.x2(cdid==mktid,:);
pcoefiM = pcoefi(cdid==mktid,:);
deltanpM = deltanp(cdid==mktid);
dfullM = vars.dfull(cdid==mktid,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing prices and shares.
% - appears robust to starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prices for coalition firms
if spec.bysize==0
    price(coalid==1) = pnashx(coalid==1) + sm;
elseif spec.bysize==1
    price(coalid==1 & sizeid==1) = pnashx(coalid==1 & sizeid==1) + sm(1);
    price(coalid==1 & sizeid==2) = pnashx(coalid==1 & sizeid==2) + sm(2);
    price(coalid==1 & sizeid==3) = pnashx(coalid==1 & sizeid==3) + sm(3);
elseif spec.bysize==2
    price(coalid==1 & sizeid==1) = pnashx(coalid==1 & sizeid==1) + sm(1);
    price(coalid==1 & sizeid==2) = pnashx(coalid==1 & sizeid==2) + sm(1);
    price(coalid==1 & sizeid==3) = pnashx(coalid==1 & sizeid==3) + sm(2);
end

% Adjustment for Modelo in coalition with separate m (counterfactual)
if spec.modsm==1
    price(brndid==5 | brndid==6) = pnashx(brndid==5 | brndid==6) + modm;
end



% Prices for fringe firms
owner   = f_ownMat(firmid(coalid==0));
p = 2*mcx(coalid==0);
popt = 1- coalid;
pfix = price(coalid==1);
f = @(p)cf_foc_partial(p,daugfile.alpha,daugfile.rho,deltanpM,mcx(coalid==0),...
    owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns);
[fpn,~,test,~] = fsolve(f,p,spec.options.cf);
price(coalid==0) = fpn;
clear owner p popt pfix fpn

% Shares
owner   = f_ownMat(firmid);
p = price;
popt = ones(size(p));
pfix = [];
[~,share,~,~] = cf_foc_partial(p,daugfile.alpha,daugfile.rho,deltanpM,mcx,...
    owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing profit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi = (price - mcx).*share.*vars.msize(cdid==mktid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtaining deviation price, shares, and profit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdev = 0;
sdev = 0;
pidev = 0;

if getdev==1
    
    [pdev,sdev,pidev] = f_getdev(mktid,price,mcx,vars,ids,daugfile,spec,restr);
    
end








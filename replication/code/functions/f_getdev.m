function [pdev,sdev,pidev] = f_getdev(mktid,price,mc,vars,ids,daugfile,spec,restr)
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the deviation prices, shares, and profit, for a
% given price vector and for a specific period*city "market" defined by
% cdid==mktid.  The output objects are matrices with one column for each
% coalition firm's deviation.
%
% Called by:
%   - f_pi_m.m
%
% Calls the following user-specified functions:
%   - f_ownMat.m
%   - cf_foc_partial.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm and coalition identity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdid = ids.cdid;
firmid = ids.firmid(cdid==mktid);
coalid = ids.coalid(cdid==mktid);
coalfirms = unique(firmid(coalid==1));

% Used in the case that derivatives are restricted
fy = unique(ids.fiscid(cdid==mktid));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demand data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x2M = vars.x2(cdid==mktid,:);

pcoefi = daugfile.pcoefi(ids.obsindemand==1,:);
deltanp = daugfile.deltanp(ids.obsindemand==1);

pcoefiM = pcoefi(cdid==mktid,:);
deltanpM = deltanp(cdid==mktid);
dfullM = vars.dfull(cdid==mktid,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtaining deviation price, shares, and profit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdev = repmat(price,[1 length(coalfirms)]);
sdev = zeros(size(pdev));
pidev = zeros(size(pdev));

countnum = 1;
for fd = coalfirms'
    
    % Option to restrict to Coors/MillerCoors to speed computation
    if restr==0 || (restr==1 && ((fy<2008 && fd==4) || (fy>2008 && fd==5)))
        
        %Deviation prices
        owner   = f_ownMat(firmid(firmid==fd));
        p = price(firmid==fd);
        popt = firmid==fd;
        pfix = price(firmid~=fd);
        f = @(p)cf_foc_partial(p,daugfile.alpha,daugfile.rho,deltanpM,mc(firmid==fd),...
            owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns);
        [fpd,~,~,~] = fsolve(f,p,spec.options.cf);
        pdev(firmid==fd,countnum) = fpd;
        clear owner p popt pfix fpn
        
        %Deviation shares
        owner   = f_ownMat(firmid);
        p = pdev(:,countnum);
        popt = ones(size(p));
        pfix = [];
        [~,sdev(:,countnum),~,~] = cf_foc_partial(p,daugfile.alpha,daugfile.rho,deltanpM,mc,...
            owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns);
        clear owner p popt pfix
        
        %Profit
        pidev(:,countnum) = (pdev(:,countnum) - mc).*sdev(:,countnum).*vars.msize(cdid==mktid);
        
    end
    
    %Iterate
    countnum = countnum+1;
    
end








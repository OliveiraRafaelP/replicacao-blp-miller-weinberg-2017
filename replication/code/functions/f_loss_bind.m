function [zero,moments,outdata,balance,bindid,gicc] = ...
    f_loss_bind(smvec,mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,modsm)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the objective function value for the imputation of
% the city-specifi supermarkups in fiscal year "fy" under the
% assumption that an ***ICC binds***.  The moments that equal zero are: 
%   (1) g(sm) = 0 , i.e., pooled slack function equals zero 
%   (2) nonlinear equations from notes
% The value of eqsolve determines whether the output "loss" is a vector of
% these moments (eqsolve==1) or a sum of squares (eqsolve==0).
% The function also returns the implied marginal costs, Bertrand prices, 
% Bertrand shares, deviation prices, and deviation shares.
%
% Called by:
%   - main_supply_bind.m
%   - f_rebalance.m
%   - results_baseanalysis.m
%   - results_cfmergers.m
%   - results_cfmmc.m
%   - results_cfanalysis.m
%   - idgraph.m
% Calls the following user-specified functions:
%   - f_impute_mc.m
%   - f_pi_m.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pnash = 0;
% snash = 0;
% mc    = 0;
% dpi   = 0;

% Supermarkups: smvec is base followed by deviations
if lspec.fullvec==1
    sm = smvec;
    sm(2:end) = smvec(2:end) + smvec(1);

% Supermarkups: smvec has base only; smdevs in smdevs
elseif lspec.fullvec==0
    sm = [smvec ; smvec+smdevs];
    
% Supermarkups: smvec has city-specific supermarkups    
elseif lspec.fullvec==2
    sm = smvec;
end


% Special supermarkup for Modelo---usually does not exist
if spec.modsm==1
    modsm = modsm;
elseif spec.modsm==0
    modsm = zeros(size(sm));
end



% Objects to help meet requirements of parfor loop (in mp function)
mktnums = 1:length(mktvec);
numProd = length(unique(ids.prodid));
temp1 = ids.coalid(ids.cdid==mktvec(1));
temp2 = ids.firmid(ids.cdid==mktvec(1));
numCoalF = length(unique(temp2(temp1==1)));
clear temp1 temp2

% Getting the fiscal year (need with lspec.icc_restrict==1)
fy = unique(ids.fiscid(ids.cdid==min(mktvec)));
    
% Shell files
mcmat = -999*ones(numProd,length(mktvec));
pnashmat = mcmat;
snashmat = mcmat;
pinashmat = mcmat;
pcollmat = mcmat;
scollmat = mcmat;
picollmat = mcmat;
pdevmat = mcmat;
sdevmat = mcmat;
pidevmat = mcmat;
giccmat = -999*ones(numCoalF,length(mktvec));
dpimat = giccmat;
dpidevmat = giccmat;

% tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  First loop is for marginal costs, Bertrand, deviation, ICCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through each of the specified CDID (city*period) within mktvec
for m=mktnums

    m;
    
    % Picking the CDID
    mktid = mktvec(m);
    
    % Matching the supermarkup
    city = unique(ids.cityid(ids.cdid==mktid));
    sm0  = sm(spec.city_in==city);
    modm0 = modsm(spec.city_in==city);
    
    % Only do calculations for cities of interest
    loc = find(spec.city_in==city);
    if lspec.skipcity(loc)==0
        
        
        % Coalagg is useful in aggregating within the market
        firmid = ids.firmid(ids.cdid==mktid);
        coalid = ids.coalid(ids.cdid==mktid);
        coalfirms = unique(firmid(coalid==1));
        temp = repmat(coalfirms',[length(firmid) 1]);
        coalagg = repmat(firmid,[1 length(coalfirms)]) == temp;
        
        % Imputation ('I') or Simulation ('S')
        if lspec.purpose=='I'
            
            % Impute marginal costs and Bertrand prices/shares
            [mc0,pnash0,snash0] = f_impute_mc(sm0,mktid,vars,ids,daugfile,spec);
            
            if spec.modsm==1 
                disp('ERROR: Imputation with Modelo Supermarkup') 
            end
            
        elseif lspec.purpose=='S'
            
            % Costs and Nash price/quantity for this market
            mc0 = idat.mc(ids.cdid==mktid);
            pnash0 = idat.pnash(ids.cdid==mktid);
            snash0 = idat.snash(ids.cdid==mktid);
                        
        end

        % Bertrand profit
        pinash0 = (pnash0 - mc0).*snash0.*vars.msize(ids.cdid==mktid);

        % PL price and shares for this market, given supermarkup sm0
        [pcoll0,scoll0,picoll0,pdev0,sdev0,pidev0] = ...
            f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,1,0,modm0);
        
        % ICC Evaluation--contribution of this CDID--only fills in
        % Coors/MillerCoors if lspec.icc_restrict==1
        temp1 = 1/(1-spec.df)*repmat(picoll0,[1 length(coalfirms)]);
        temp2 = spec.df/(1-spec.df)*repmat(pinash0,[1 length(coalfirms)]);
        temp3 = temp1 - (pidev0 + temp2);
        gicc0 = sum(temp3.*coalagg);
       
        % Storing marginal cost in shell file
        mcx = [mc0; -998*ones(numProd-length(mc0),1)];
        mcmat(:,m) = mcx;
        
        % Storing Bertrand prics and shares in shell file
        pnashx = [pnash0; -998*ones(numProd-length(mc0),1)];
        snashx = [snash0; -998*ones(numProd-length(mc0),1)];
        pinashx = [pinash0; -998*ones(numProd-length(mc0),1)];
        pnashmat(:,m) = pnashx;
        snashmat(:,m) = snashx;
        pinashmat(:,m) = pinashx;
        
        % Storing PL prices and shares in shell file
        pcollx = [pcoll0; -998*ones(numProd-length(mc0),1)];
        scollx = [scoll0; -998*ones(numProd-length(mc0),1)];
        picollx = [picoll0; -998*ones(numProd-length(mc0),1)];
        pcollmat(:,m) = pcollx;
        scollmat(:,m) = scollx;
        picollmat(:,m) = picollx;

        % Storing deviation prices and shares in shell file
        pdevx = [min(pdev0,[],2); -998*ones(numProd-length(mc0),1)];
        sdevx = [max(sdev0,[],2); -998*ones(numProd-length(mc0),1)];
        pidevx = [max(pidev0,[],2); -998*ones(numProd-length(mc0),1)];
        pdevmat(:,m) = pdevx;
        sdevmat(:,m) = sdevx;
        pidevmat(:,m) = pidevx;
        
        % Storing ICC contribution in shell file
        giccmat(:,m) = gicc0';
        
    elseif lspec.skipcity(loc)==1
        
        filler = zeros(size(ids.cdid(ids.cdid==mktid)));
        fillerx = [filler; -998*ones(numProd-length(filler),1)];
        mcmat(:,m) = fillerx;
        pnashmat(:,m) = fillerx;
        snashmat(:,m) = fillerx;
        pcollmat(:,m) = fillerx;
        scollmat(:,m) = fillerx;
        pdevmat(:,m) = fillerx;
        sdevmat(:,m) = fillerx;        
        picollmat(:,m) = fillerx;
        pinashmat(:,m) = fillerx;
        pidevmat(:,m) = fillerx;
        %pidevmat(:,m,:) = repmat(fillerx,[1 numCoalF]);
        
    end
    
end
    
% Marginal cost vector
mc = mcmat(:);
mc = mc(mc~= -998);

% Bertrand prices and shares and profit
pnash = pnashmat(:);
snash = snashmat(:);
pinash = pinashmat(:);
pnash = pnash(pnash~= -998);
snash = snash(snash~= -998);
pinash = pinash(pinash~= -998);

% PL prices and shares and profit
pcoll = pcollmat(:);
scoll = scollmat(:);
picoll = picollmat(:);
pcoll = pcoll(pcoll~= -998);
scoll = scoll(scoll~= -998);
picoll = picoll(picoll~= -998);

% Deviation prices and shares and profit
pdev = pdevmat(:);
sdev = sdevmat(:);
pidev = pidevmat(:);
pdev = pdev(pdev~= -998);
sdev = sdev(sdev~= -998);
pidev = pidev(pidev~= -998);

% Output structure
outdata.mc=mc;
outdata.pnash = pnash;
outdata.snash = snash;
outdata.pcoll = pcoll;
outdata.scoll = scoll;
outdata.pdev  = pdev;
outdata.sdev  = sdev;
outdata.pinash = pinash;
outdata.picoll = picoll;
outdata.pidev = pidev;

%toc/60

% Slack values, "most violating" slack value, identity of binder
gicc = sum(giccmat,2);
if lspec.icc_restrict==0
    
    gbind = min(gicc);
    binder = gicc == gbind;
    
elseif lspec.icc_restrict==1
    
    if fy<2008
        binder = [0;1;0];
    elseif fy>2008
        binder = [0;1];
    end
    gbind = gicc(binder==1);
    
elseif lspec.icc_restrict==2
    
    if fy<2008
        binder = [0;1;0;0];
    elseif fy>2008
        binder = [0;1;0];
    end
    gbind = gicc(binder==1);
    
elseif lspec.icc_restrict==3

    binder = [0;1];
    gbind = gicc(binder==1);
    
end

% Tie-break rule which should not be necessary but just in case
if sum(binder)>1
   for i=1:length(binder)
       if sum(binder)>1
           binder(i)=0;
       end
   end
end

% A bit of hard-coding here: firmid of binding firm
if find(binder==1)==1
    bindid=1;
elseif find(binder==1)==2 && fy < 2008
    bindid=4;
else
    bindid=5;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Second loop for numerical derivatives, then balance evaluations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Only if requested!
if lspec.getder==1
    
    % Need for slotting in data
    cdid2 = ids.cdid(ismember(ids.cdid,mktvec)==1);
    
    % Loop through each of the specified CDID (city*period) within mktvec
    for m=mktnums
        
        % Picking the CDID
        mktid = mktvec(m);
        
        % Matching the supermarkup
        city = unique(ids.cityid(ids.cdid==mktid));
        sm0  = sm(spec.city_in==city);
        modm0 = modsm(spec.city_in==city);
            
        % Only do calculations for cities of interest
        loc = find(spec.city_in==city);
        if lspec.skipcity(loc)==0
            
            % Grabbing relevant data
            pnash0 = pnash(cdid2==mktid);
            mc0 = mc(cdid2==mktid);
            
            % Coalagg is useful in aggregating within the market
            firmid = ids.firmid(ids.cdid==mktid);
            coalid = ids.coalid(ids.cdid==mktid);
            coalfirms = unique(firmid(coalid==1));
            temp = repmat(coalfirms',[length(firmid) 1]);
            coalagg = repmat(firmid,[1 length(coalfirms)]) == temp;
            
            % Perturbing the the supermarkup
            eps = 0.01;
            sm1 = sm0 + eps;
            [~,~,pi_b,~,~,pidev_b]=f_pi_m(sm1,mktid,pnash0,mc0,vars,ids,daugfile,spec,1,lspec.icc_restrict,modm0);
            
            sm2 = sm0 - eps;
            [~,~,pi_a,~,~,pidev_a]=f_pi_m(sm2,mktid,pnash0,mc0,vars,ids,daugfile,spec,1,lspec.icc_restrict,modm0);
            
            % Calculating two-sided numerical derivative
            temp1 = repmat(pi_b - pi_a,[1 length(coalfirms)]);
            temp2 = pidev_b - pidev_a;
            dpi = sum(temp1.*coalagg);
            dpidev = sum(temp2.*coalagg);
            
            % Storing derivative of PLE profit in shell file
            dpimat(:,m) = dpi';
            
            % Storing derivative of deviation profit in shell file
            dpidevmat(:,m) = dpidev';
            
        end
        
    end
    
end

%toc

% For the counterfactual with MillerCoors as the leader
if lspec.mcleader==1
    nn = 2;
else
    nn = 1;
end


% Market-specific calculations, if requested
temp1 = ids.cityid(ids.cdindex);
citymatch = temp1(mktvec);  

if lspec.getder==1 && length(unique(citymatch))>1
    leader_numer = zeros(length(unique(citymatch)),1);
    binder_denom = zeros(length(unique(citymatch)),1);
    for r = 1:length(leader_numer)
        leader_numer(r) = sum(dpimat(nn,citymatch==spec.city_in(r)));
        binder_denom(r) = sum(dpidevmat(binder==1,citymatch==spec.city_in(r)))  ...
                        - sum(dpimat(binder==1,citymatch==spec.city_in(r)));
    end
    bkeep = leader_numer ./ binder_denom;
    balance.rat = bkeep ./ bkeep(1);
    balance.lev = bkeep - bkeep(1);
    balance.num = leader_numer;
    balance.den = binder_denom;
    
    % The following four lines are important for the skipcity flag to work
    balance.rat(lspec.skipcity==1)=1;
    balance.lev(lspec.skipcity==1)=0;
    balance.num(lspec.skipcity==1)=0;
    balance.den(lspec.skipcity==1)=0;
    
elseif lspec.getder==0 || length(unique(citymatch))==1
    balance = 0;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Third loop for numerical derivatives of special Modelo supermarkup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only if requested!
if lspec.getder==1 && spec.modsm==1
    
    % Need for slotting in data
    cdid2 = ids.cdid(ismember(ids.cdid,mktvec)==1);
    
    % Loop through each of the specified CDID (city*period) within mktvec
    for m=mktnums
        
        % Picking the CDID
        mktid = mktvec(m);
        
        % Matching the supermarkup
        city = unique(ids.cityid(ids.cdid==mktid));
        sm0  = sm(spec.city_in==city);
        modm0 = modsm(spec.city_in==city);
            
        % Only do calculations for cities of interest
        loc = find(spec.city_in==city);
        if lspec.skipcity(loc)==0
            
            % Grabbing relevant data
            pnash0 = pnash(cdid2==mktid);
            mc0 = mc(cdid2==mktid);
            
            % Coalagg is useful in aggregating within the market
            firmid = ids.firmid(ids.cdid==mktid);
            coalid = ids.coalid(ids.cdid==mktid);
            coalfirms = unique(firmid(coalid==1));
            temp = repmat(coalfirms',[length(firmid) 1]);
            coalagg = repmat(firmid,[1 length(coalfirms)]) == temp;
            
            % Perturbing the supermarkup
            eps = 0.01;
            modm1 = modm0 + eps;
            [~,~,pi_b,~,~,pidev_b]=f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,1,lspec.icc_restrict,modm1);
            
            modm2 = modm0 - eps;
            [~,~,pi_a,~,~,pidev_a]=f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,1,lspec.icc_restrict,modm2);
            
            % Calculating two-sided numerical derivative
            temp1 = repmat(pi_b - pi_a,[1 length(coalfirms)]);
            temp2 = pidev_b - pidev_a;
            dpi = sum(temp1.*coalagg);
            dpidev = sum(temp2.*coalagg);
            
            % Storing derivative of PLE profit in shell file
            dpimat(:,m) = dpi';
            
            % Storing derivative of deviation profit in shell file
            dpidevmat(:,m) = dpidev';
            
        end
        
    end
    
end

%toc


% Market-specific calculations, if requested
temp1 = ids.cityid(ids.cdindex);
citymatch = temp1(mktvec);

if lspec.getder==1 && length(unique(citymatch))>1 && spec.modsm==1
    leader_numer = zeros(length(unique(citymatch)),1);
    binder_denom = zeros(length(unique(citymatch)),1);
    for r = 1:length(leader_numer)
        leader_numer(r) = sum(dpimat(1,citymatch==spec.city_in(r)));
        binder_denom(r) = sum(dpidevmat(binder==1,citymatch==spec.city_in(r)))  ...
                        - sum(dpimat(binder==1,citymatch==spec.city_in(r)));
    end
    balance.ratmod = (leader_numer ./ binder_denom) ./ bkeep(1);
    balance.levmod = (leader_numer ./ binder_denom) - bkeep(1);
    balance.nummod = leader_numer;
    balance.denmod = binder_denom;
    
    % The following four lines are important for the skipcity flag to work
    balance.ratmod(lspec.skipcity==1)=1;
    balance.levmod(lspec.skipcity==1)=0;
    balance.nummod(lspec.skipcity==1)=0;
    balance.denmod(lspec.skipcity==1)=0;
    
end

%balance.num';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Moments and Loss 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moments that equal zero for PLE supermarkup
if lspec.getder==1 && length(unique(citymatch))>1
    moments = [gbind/100 ; lspec.scale*balance.lev];
elseif lspec.getder==0 || length(unique(citymatch))==1
    moments = gbind;
end

% Loss function 
zero = sum(moments.^2);



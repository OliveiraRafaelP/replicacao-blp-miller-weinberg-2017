function [zero,gbind,mc,pnash,snash,balance,bindid,gicc] = ...
    f_loss_bind_x(smbase,mktvec,vars,ids,idat,daugfile,spec,lspecx,smdevs)


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
%   - main_supply_bind_x.m
%   - f_rebalance_x.m
% Calls the following user-specified functions:
%   - f_impute_mc.m
%   - f_pi_m.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pnash = 0;
% snash = 0;
% mc    = 0;
% dpi   = 0;

% lspecx = lspec_2;

% Vector/matrix of supermarkups
sm = smbase+smdevs;

% Objects to help meet requirements of parfor loop (in mp function)
mktnums = 1:length(mktvec);
numProd = length(unique(ids.prodid));
temp1 = ids.coalid(ids.cdid==mktvec(1));
temp2 = ids.firmid(ids.cdid==mktvec(1));
numCoalF = length(unique(temp2(temp1==1)));
clear temp1 temp2
 
% Shell files
mcmat = zeros(numProd,length(mktvec));
pnashmat = mcmat;
snashmat = mcmat;
pimat = mcmat;
pidevmat = repmat(mcmat,[1 1 numCoalF]);
giccmat = zeros(numCoalF,length(mktvec));
dpimat = repmat(giccmat,[1 1 3]);
dpidevmat = repmat(giccmat,[1 1 3]);

% tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  First loop is for marginal costs, Bertrand, deviation, ICCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through each of the specified CDID (city*period) within mktvec
for m=mktnums

    % Picking the CDID
    mktid = mktvec(m);
     
    % Matching the supermarkup
    city = unique(ids.cityid(ids.cdid==mktid));
    if spec.pool==0
        sm0 = sm;
    elseif spec.pool==1
        if spec.bysize==0
            sm0 = sm(spec.city_in==city);
        else
            sm0 = sm(spec.city_in==city,:);
        end
    end
    
    % Only do calculations for cities of interest
    loc = find(spec.city_in==city);
    if lspecx.skipcity(loc)==0
        
        % Coalagg is useful in aggregating within the market
        firmid = ids.firmid(ids.cdid==mktid);
        coalid = ids.coalid(ids.cdid==mktid);
        coalfirms = unique(firmid(coalid==1));
        temp = repmat(coalfirms',[length(firmid) 1]);
        coalagg = repmat(firmid,[1 length(coalfirms)]) == temp;
        
        
        % Imputation ('I') or Simulation ('S')
        if lspecx.purpose=='I'
            
            % Impute marginal costs and Bertrand prices/shares
            [mc0,pnash0,snash0] = f_impute_mc(sm0,mktid,vars,ids,daugfile,spec);
            
        elseif lspecx.purpose=='S'
            
            % Costs and Nash price/quantity for this market
            mc0 = idat.mc(ids.cdid==mktid);
            pnash0 = idat.pnash(ids.cdid==mktid);
            snash0 = idat.snash(ids.cdid==mktid);
            
        end
        
        % Bertrand profit
        pinash0 = (pnash0 - mc0).*snash0.*vars.msize(ids.cdid==mktid);
        
        % PL price and shares for this market, given supermarkup sm0
        [pcoll0,scoll0,picoll0,pdev0,sdev0,pidev0] = ...
            f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,1,0);
        
        

        % ICC Evaluation--contribution of this CDID--only fills in
        % Coors/MillerCoors if lspec.icc_restrict==1
        temp1 = 1/(1-spec.df)*repmat(picoll0,[1 length(coalfirms)]);
        temp2 = spec.df/(1-spec.df)*repmat(pinash0,[1 length(coalfirms)]);
        temp3 = temp1 - (pidev0 + temp2);
        gicc0 = sum(temp3.*coalagg);
        
        % Careful w/ icc_restrict==1 (!!)
        if lspecx.icc_restrict==1
            if spec.fy<2008
                gicc0([1;3]) = [0;0];
            elseif spec.fy>2008
                gicc0(1) = 0;
            end
        end
        
        % Storing marginal cost in shell file
        mcx = [mc0; -998*ones(numProd-length(mc0),1)];
        mcmat(:,m) = mcx;
        
        % Storing Bertrand prics and shares in shell file
        pnashx = [pnash0; -998*ones(numProd-length(mc0),1)];
        snashx = [snash0; -998*ones(numProd-length(mc0),1)];
        pnashmat(:,m) = pnashx;
        snashmat(:,m) = snashx;
        
        % Storing profits for use with numberical differences
        pi0x = [picoll0; -998*ones(numProd-length(mc0),1)];
        pidev0x = [pidev0; -998*ones(numProd-length(mc0),numCoalF)];
        pimat(:,m) = pi0x;
        pidevmat(:,m,:) = pidev0x;
       
        % Storing ICC contribution in shell file
        giccmat(:,m) = gicc0';
        
    elseif lspecx.skipcity(loc)==1
        
        filler = zeros(size(ids.cdid(ids.cdid==mktid)));
        fillerx = [filler; -998*ones(numProd-length(filler),1)];
        mcmat(:,m) = fillerx;
        pnashmat(:,m) = fillerx;
        snashmat(:,m) = fillerx;
        pimat(:,m) = fillerx;
        pidevmat(:,m,:) = repmat(fillerx,[1 numCoalF]);
        
    end
    
end

% Marginal cost vector
mc = mcmat(:);
mc = mc(mc~= -998);

% Bertrand prices and shares
pnash = pnashmat(:);
snash = snashmat(:);
pnash = pnash(pnash~= -998);
snash = snash(snash~= -998);

% PL and Binder's Deviation profit
picoll = pimat(:);
picoll = picoll(picoll~= -998);
pidev = reshape(pidevmat,[],size(pidevmat,3),1);
temp = pidev(:,1)==-998;
pidev(temp==1,:) = [];

% Slack values, "most violating" slack value, identity of binder
gicc = sum(giccmat,2);
if lspecx.icc_restrict==0
    gbind = min(gicc);
    binder = gicc == gbind;
elseif lspecx.icc_restrict==1
    if spec.fy<2008
        binder = [0;1;0];
    elseif spec.fy>2008
        binder = [0;1];
    end
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
elseif find(binder==1)==2 && spec.fy < 2008
    bindid=4;
else
    bindid=5;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Second loop for numerical derivatives, then balance evaluations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Foor looping
if spec.bysize==0
    maxloop = 1;
elseif spec.bysize==1
    maxloop = 3;
elseif spec.bysize==2
    maxloop = 2;
end


% Only if requested!
if lspecx.getder==1
    
    % Need for slotting in data
    cdid2 = ids.cdid(ismember(ids.cdid,mktvec)==1);

    
    % Loop through each of the specified CDID (city*period) within mktvec
    for m=mktnums
        
        % Picking the CDID
        mktid = mktvec(m);
        
        % Matching the supermarkup
        city = unique(ids.cityid(ids.cdid==mktid));
        if spec.pool==0
            sm0 = sm;
        elseif spec.pool==1
            if spec.bysize==0
                sm0 = sm(spec.city_in==city);
            else
                sm0 = sm(spec.city_in==city,:);
            end
        end
        
        % Only do calculations for cities of interest
        loc = find(spec.city_in==city);
        if lspecx.skipcity(loc)==0
            
            % This allows us to get the derivative WRT each sm
            for sid = 1:maxloop
                
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
                eps = 1e-3;  % 0.001;
                sm1 = sm0;
                sm1(sid) = sm1(sid) + eps;
                [~,~,pi_b,~,~,pidev_b]=f_pi_m(sm1,mktid,pnash0,mc0,vars,ids,daugfile,spec,1,lspecx.icc_restrict);
                
                % Using single or double-sided numerical differentiation
                if lspecx.nd2==1
                    sm2 = sm0;
                    sm2(sid) = sm2(sid) - eps;
                    [~,~,pi_a,~,~,pidev_a]=f_pi_m(sm2,mktid,pnash0,mc0,vars,ids,daugfile,spec,1,lspecx.icc_restrict);
                    gap = 1;
                elseif lspecx.nd2==0
                    
                    pi_a = picoll(cdid2==mktid);
                    pidev_a = pidev(cdid2==mktid,:);
                    gap = 2;
                    
                end
                
                % Calculating two-sided numerical derivative
                temp1 = repmat(pi_b - pi_a,[1 length(coalfirms)]);
                temp2 = pidev_b - pidev_a;
                dpi = gap*sum(temp1.*coalagg);
                dpidev = gap*sum(temp2.*coalagg);
                
                % Storing derivative of PLE profit in shell file
                dpimat(:,m,sid) = dpi';
                
                % Storing derivative of deviation profit in shell file
                dpidevmat(:,m,sid) = dpidev';
                
            end
        end
    end
end



% Market-specific calculations, if requested
if lspecx.getder==1 
    
    temp1 = ids.cityid(ids.cdindex);
    citymatch = temp1(mktvec);

    % Shell files
    leader_numer = ones(length(unique(citymatch)),3);
    binder_denom = ones(length(unique(citymatch)),3);
    
    % Filling in shell files---careful here with indexing
    for sid = 1:maxloop
        if spec.pool==1
            for r = 1:size(leader_numer,1)
                leader_numer(r,sid) = sum(dpimat(1,citymatch==spec.city_in(r),sid));
                binder_denom(r,sid) = sum(dpidevmat(binder==1,citymatch==spec.city_in(r),sid) ...
                    - dpimat(binder==1,citymatch==spec.city_in(r),sid) );
            end
        elseif spec.pool==0
            leader_numer(sid) = sum(dpimat(1,:,sid));
            binder_denom(sid) = sum(dpidevmat(binder==1,:,sid)-dpimat(binder==1,:,sid) );
        end
    end
    
    % Eliminating extra dimensions
    if spec.bysize==0
        leader_numer(:,2:3) = [];
        binder_denom(:,2:3) = [];
    elseif spec.bysize==2
        leader_numer(:,3) = [];
        binder_denom(:,3) = [];
    end
    
    % Creating balance terms
    btemp = leader_numer ./ binder_denom;
    balance.rat = btemp ./ btemp(1);
    balance.lev = btemp - btemp(1);
    balance.num = leader_numer;
    balance.den = binder_denom;
    
    % The following two lines are important for the skipcity flag to work
    balance.rat(isnan(balance.rat))=1;
    balance.lev(isnan(balance.lev))=0;    
    balance.num(lspecx.skipcity==1,:)=0;
    balance.den(lspecx.skipcity==1,:)=0;    
    
elseif lspecx.getder==0 
    
    balance = 0;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining loss function materials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loss function ; poentially modified below
zero = gbind^2;



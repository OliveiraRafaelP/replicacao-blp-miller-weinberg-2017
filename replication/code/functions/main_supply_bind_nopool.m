function main_supply_bind_nopool(path,spec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function imputes supermarkups and marginal costs under assumption
% of **binding** ICCs.  The input arguments specify a discount factor
% ("df"). The model features region-specific supermarkups with ICCs that
% are *not* pooled across cities. Thus, the identity of the binding firm
% can (and will) differ depending on the city.
%
% Called by:
%   - main_v7.m
% Calls the following user-specified functions:
%   - main_data.m
%   - f_loss_bind.m
%   - f_rebalance_sp.m
% Saved results are loaded in:
%   - results_costregs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data required for imputation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Loading demand results
cd(strcat(path.data2))
load daugfile daugfile
cd(strcat(path.code1))

% Main data
[vars,ids] = main_data(path,spec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shell files
mc = 0*repmat(vars.p_jt,[1 4]);
pnash = mc;
snash = mc;
sm_bind = zeros(4,39);
fval = zeros(4,39);
bindid = zeros(4,39);

% Optimization options: only use tight
optsimprb_2 = optimset('GradObj','off','MaxIter',50000,'MaxFunEvals',100000,...
    'Display','off','TolFun',1e-4,'TolX',1e-4);

% Initial imputation options
lspec.fullvec=0;
lspec.icc_restrict=0;
lspec.getder=0;
lspec.scale = 100;      % won't matter but need something
lspec.purpose='I';
lspec.skipcity = zeros(length(spec.city_in),1);
lspec.mcleader = 0;

% Necessary to run (I think)
smdevs = zeros(length(spec.city_in)-1,1);

% Starting values
startval = [0.75 0.75 1.25 1.25];
startval = startval + spec.df*startval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system of equations using the obtained starting value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For slotting in results across fiscal years;  starting timer
counter=1;
tic

% Looping over the fiscal years with complete data
for fy = [2006,2007,2010,2011]
    
    disp(fy);
    
    numProd = length(unique(ids.prodid));
    numMon = max(ids.montid);
    mcmat = -999*ones(numProd*numMon,39);
    pnashmat = mcmat;
    snashmat = mcmat;
    
    %Run city-by-city; use mp across cities (highest level)
    for c = spec.city_in'
     
        disp(c);
        
        % Identify corresponding CDID codes
        mktvec =  unique(ids.cdid(ids.fiscid==fy & ids.cityid==c));
        
        % Imputing constrained supermarkup--simplex due to multiple mins
        fb = @(sm)f_loss_bind(sm,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs,[]);
        [smb0,fval0] = fminsearch(fb,startval(counter),optsimprb_2);
        
        % Starting values can cause convergence to zero
        if smb0 < 0.01
            disp('Convergence to zero: scaling down start');
            [smb0,fval0] = fminsearch(fb,0.50*startval(counter),optsimprb_2);
        end
        if smb0 < 0.01
            disp('Convergence to zero: scaling down start again');
            [smb0,fval0] = fminsearch(fb,0.25*startval(counter),optsimprb_2);
        end
        if smb0 < 0.01
            disp('Convergence to zero: scaling up start');
            [smb0,fval0] = fminsearch(fb,1.50*startval(counter),optsimprb_2);
        end
        if smb0 < 0.01
            disp('Final Convergence to zero. Check:')
            disp([spec.df*100 fy c]);
        end
        
       
        % Imputing other statistics
        [~,~,outdata,~,bind0,~] = f_loss_bind(smb0,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs,[]);
        if fy>2008 && bind0==4
            bind0=5;
        end
        
        %Slotting in results
        mcx = [outdata.mc; -998*ones(numMon*numProd-length(outdata.mc),1)];
        pnashx = [outdata.pnash; -998*ones(numMon*numProd-length(outdata.mc),1)];
        snashx = [outdata.snash; -998*ones(numMon*numProd-length(outdata.mc),1)];
        mcmat(:,c) = mcx;
        pnashmat(:,c) = pnashx;
        snashmat(:,c) = snashx;
        sm_bind(counter, c) = smb0;
        fval(counter, c) = fval0;
        bindid(counter, c) = bind0;
        
        %disp([c smb0 fval0])
        
    end
    
    %Extracting results matrices
    for c=spec.city_in'
        temp = mcmat(:,c);
        temp = temp(temp~=-998);
        mc(ids.cityid==c & ids.fiscid==fy) = temp;
        
        temp = pnashmat(:,c);
        temp = temp(temp~=-998);
        pnash(ids.cityid==c & ids.fiscid==fy) = temp;
        
        temp = snashmat(:,c);
        temp = temp(temp~=-998);
        snash(ids.cityid==c & ids.fiscid==fy) = temp;
    end
       
    %Updating counter & moving to next fiscal year
    counter = counter+1;
    
end

shours = toc/60/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transposing results and keeping the relevant cities and fiscal years
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sm = sm_bind(:, spec.city_in)';
fval = fval(:, spec.city_in)';
bindid = bindid(:, spec.city_in)';

%%%%%%%%%%%%%%%%%%
% Storing results
%%%%%%%%%%%%%%%%%%

cd(strcat(path.data2,'/',spec.ssubfolder))
    
    save sres_bind_nopool mc pnash snash sm shours bindid fval

cd(strcat(path.code1))


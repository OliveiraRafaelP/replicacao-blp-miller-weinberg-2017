function results_cfmmc(path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes equilibrium in the counterfactual featuring no
% multimarket contact.
%
% Called by:
%   - main_v6.m
% Calls the following user-specified functions:
%   - main_data.m
%   - main_spec.m
%   - f_loss_bind.m
% Loads results created in the functions:
%   - combine_imputed.m
% Saved results are loaded in:
%   - results_cfanalysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Baseline specification
dfx = 26;
df = dfx/100;
spec = main_spec(df);
spec.bysize = 0;  

% Grabbing data and defining fiscnum
[vars,ids] = main_data(path,spec);

% Demand results
cd(strcat(path.data2))
load daugfile daugfile;
cd(path.code1)


% Results at selected discount factor
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc pnash snash sm_06 sm_07 sm_10 sm_11 ;
cd(path.code1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input data
idat.mc = mc;
idat.pnash = pnash;
idat.snash = snash;
clear mc pnash snash

% Simulation options
lspec.getder = 0;       
lspec.fullvec = 0;      
lspec.icc_restrict = 0; % restrict ICC analysis to Coors/MillerCoors
lspec.scale = 100;      
lspec.purpose = 'S';    
lspec.skipcity = zeros(length(spec.city_in),1);
lspec.mcleader = 0;

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

% Necessary to run (I think)
smdevs = zeros(length(spec.city_in)-1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system of equations using the obtained starting value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For slotting in results across fiscal years;  starting timer
counter=1;
tic

% Looping over the fiscal years with complete data
for fy = [2006,2007,2010,2011]
    
    % Starting values
    if fy==2006
        startval = sm_06;
    elseif fy==2007
        startval = sm_07;
    elseif fy==2010
        startval = sm_10;
    elseif fy==2011
        startval = sm_11;
    end
    
    disp(fy);
    
    numProd = length(unique(ids.prodid));
    numMon = max(ids.montid);
    mcmat = -999*ones(numProd*numMon,39);
    pnashmat = mcmat;
    snashmat = mcmat;
    
    %Run city-by-city; use mp across cities (highest level)
    for c = spec.city_in'
     
        disp(c);
        
        % Identify corresponding CDID codes, and starting value
        mktvec =  unique(ids.cdid(ids.fiscid==fy & ids.cityid==c));
        sv = startval(spec.city_in==c);
        
        % Imputing constrained supermarkup--simplex 
        fb = @(sm)f_loss_bind(sm,mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,[]);
        [smb0,fval0] = fminsearch(fb,sv,optsimprb_2);
        

        % Starting values can cause convergence to zero
        if smb0 < 0.01
            disp('Convergence to zero: scaling down start');
            [smb0,fval0] = fminsearch(fb,0.50*sv,optsimprb_2);
        end
        if smb0 < 0.01
            disp('Convergence to zero: scaling down start again');
            [smb0,fval0] = fminsearch(fb,0.25*sv,optsimprb_2);
        end
        if smb0 < 0.01
            disp('Convergence to zero: scaling up start');
            [smb0,fval0] = fminsearch(fb,sv,optsimprb_2);
        end
        if smb0 < 0.01
            disp('Final Convergence to zero. Check:')
            disp([spec.df*100 fy c]);
        end
        
       
        % Imputing other statistics
        [~,~,outdata,~,bind0,~] = f_loss_bind(smb0,mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,[]);
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

cd(strcat(path.data2))
    
    save cfscen_nommc mc pnash snash sm shours bindid fval

cd(strcat(path.code1))


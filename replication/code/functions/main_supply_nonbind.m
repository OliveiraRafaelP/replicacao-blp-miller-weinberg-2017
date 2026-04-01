function main_supply_nonbind(path,spec)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function imputes supermarkups and marginal costs under assumption
% of non-binding ICCs.  The four fiscal years for which we have full data
% are considered: 2006, 2007, 2010, 2011. 
%
% Called by:
%   - main_v7.m
% Calls the following user-specified functions:
%   - main_data.m
%   - f_loss_nonbind.m
% Saved results are loaded in:
%   - results_costregs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data required for imputation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Loading demand results
cd(strcat(path.data2))
load daugfile
cd(strcat(path.code1))
    
% Main data 
[vars,ids] = main_data(path,spec);


%%%%%%%%%%%%%%%%%%
% Non-binding ICC
%%%%%%%%%%%%%%%%%%

tic

% Matrix of fiscal year * regions to consider
fyrmat = zeros(length(spec.city_in),4);
fyrmat(:,1) = unique(ids.fisccity(ids.fiscid==2006));
fyrmat(:,2) = unique(ids.fisccity(ids.fiscid==2007));
fyrmat(:,3) = unique(ids.fisccity(ids.fiscid==2010));
fyrmat(:,4) = unique(ids.fisccity(ids.fiscid==2011));

% Shell files
mc = 0*repmat(vars.p_jt,[1 4]);
pnash = mc;
snash = mc;
sm_unc = zeros(size(fyrmat));
fval_unc = sm_unc;
start = 3;

if spec.bysize==2
   sm_unc = repmat(sm_unc,[1 1 2]);
   start = [start;start];
end

% Looping over the fiscal years with complete data
for fy = 1:4
   
     fyrlist = fyrmat(:,fy);
    
     % Looping over the CDID = regions*periods individually
     counter=1;
     for fyr = unique(fyrlist)'
         
         %Imputing unconstrained supermarkup--BFGS much faster than simplex
         fu = @(sm)f_loss_nonbind(sm,fyr,vars,ids,daugfile,spec,1);
         [smu0,fval0] = fminunc(fu,start);

         
         %Imputing marginal costs (and Bertrand price/shares)
         [~,mc0,pnash0,snash0] = f_loss_nonbind(smu0,fyr,vars,ids,daugfile,spec,0);
         
         %Slotting in results
         mc(ids.fisccity==fyr,fy) = mc0;
         pnash(ids.fisccity==fyr,fy) = pnash0;
         snash(ids.fisccity==fyr,fy) = snash0;
         sm_unc(counter,fy,:) = smu0;
         fval_unc(counter,fy) = fval0;
         
         counter = counter+1;
         
         fyr
         
     end
    
end


caltime_unc = toc / 60 / 60;


mc = sum(mc,2);
pnash = sum(pnash,2);
snash = sum(snash,2);


%%%%%%%%%%%%%%%%%%
% Storing results
%%%%%%%%%%%%%%%%%%

cd(strcat(path.data2))
%mkdir(spec.sfolder)

if spec.bysize==0
    
    save sres_nonbind mc pnash snash caltime_unc sm_unc fval_unc
    
elseif spec.bysize==2
    
    save sres_nonbind_x mc pnash snash caltime_unc sm_unc fval_unc
    
end

cd(strcat(path.code1))


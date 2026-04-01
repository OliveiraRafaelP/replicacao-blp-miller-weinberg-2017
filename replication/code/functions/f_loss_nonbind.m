function [loss,mc,pnash,snash] = ...
    f_loss_nonbind(sm,fyc,vars,ids,daugfile,spec,getder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the objective function value for the imputation of
% the supermarkup in fiscal year * city combination "fyc" under the
% assumption that an ***ICC does not bind***. The objective function value 
% to be minimized is given by the square of d pi / d m.  The function also
% returns the implied marginal costs, Bertrand prices, and Bertrand shares.
%
% Called by:
%   - main_supply_nonbind.m
% Calls the following user-specified functions:
%   - f_impute_mc.m
%   - f_pi_m.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pnash = 0;
snash = 0;
mc    = 0;
dpi   = 0;
if spec.bysize==2
    dpi = [dpi;dpi];
end

% Loop through each of the four CDID city*period containted within fyc
for mktid = unique(ids.cdid(ids.fisccity==fyc))'
    
   %Impute marginal costs and profit
   [mc0,pnash0,snash0] = f_impute_mc(sm,mktid,vars,ids,daugfile,spec);
   
   %Will delete first element later
   mc    = [mc ; mc0];
   pnash = [pnash ; pnash0];
   snash = [snash ; snash0];

   %For numerical derivatives
   eps = 0.01;   
   leader = ids.leadid(ids.cdid==mktid);
   
   %Numerically evaluating dpi / dm if requested
   if getder==1
       
       if spec.bysize==0
           
           sm0 = sm - eps;
           [~,~,pi_a,~,~] = f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,0,1);
           
           sm0 = sm + eps;
           [~,~,pi_b,~,~] = f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,0,1);
           
           pi_change = pi_b - pi_a;
           dpi = dpi + sum( pi_change(leader==1));
           
       elseif spec.bysize==2
           
           sm0 = sm;
           sm0(1) = sm0(1) - eps;
           [~,~,pi_a,~,~] = f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,0,1);
           
           sm0(1) = sm0(1) + 2*eps;
           [~,~,pi_b,~,~] = f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,0,1);
           
           pi_change = pi_b - pi_a;
           dpi(1) = dpi(1) + sum( pi_change(leader==1));
           
           sm0 = sm;
           sm0(2) = sm0(2) - eps;
           [~,~,pi_a,~,~] = f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,0,1);
           
           sm0(2) = sm0(2) + 2*eps;
           [~,~,pi_b,~,~] = f_pi_m(sm0,mktid,pnash0,mc0,vars,ids,daugfile,spec,0,1);
           
           pi_change = pi_b - pi_a;
           dpi(2) = dpi(2) + sum( pi_change(leader==1));
           
       end
   end
   
end

%Output files
mc = mc(2:end);
pnash = pnash(2:end);
snash = snash(2:end);
loss = sum(dpi.^2);


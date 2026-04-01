function [smout,balz2,modsmout] = ...
    f_rebalance(smvec0,mktvec,vars,ids,idat,daugfile,spec,lspec,rbspec,modsm0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the region-specific supermarkup deviations that
% satify the FOCs, given a base supermarkup. 
%
% Called by:
%   - main_supply_bind.m
%   - results.cfmergers.m
% Calls the following user-specified functions:
%   - f_loss_bind.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lspec_2 = lspec;
lspec_2.fullvec = 1;
lspec_2.getder = 1;
    
[~,~,~,bal0,~,~] = f_loss_bind(smvec0,mktvec,vars,ids,idat,daugfile,spec,lspec_2,[],modsm0);


% smvec,mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,modsm

    %smvec,mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,modsm)
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebalancing Procedure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


baly = bal0.rat;
badval = bal0.num(2:end) <0 ;
smbase0 = smvec0(1);
smdevsz = smvec0(2:end);
disp(baly');

modsmz = modsm0;
if spec.modsm==1
    modbaly = bal0.ratmod;
    disp(modbaly')
end

tc = 0;
    


while tc < 1000    

    % Evaluate state
    smdevsz2 = smdevsz;
    z = abs(baly-1); 
    z2 = z(2:end);
    w = baly-1;
    w(w>0.5)=0.5;
    w(w<-0.5)=-0.5;
    w2 = w(2:end);
    
    if spec.modsm==1
        
        modsmz2 = modsmz;
        modz = abs(modbaly-1);
        modw = modbaly'-1;
        modw(modw>0.5)=0.5;
        modw(modw<-0.5)=-0.5;
        
    end
    
    
    % Adjust city-skipper
    if spec.pool==1
        if spec.modsm==0
            lspec_2.skipcity = z < rbspec.thresh;
        elseif spec.modsm==1 && modz(1) < rbspec.thresh
            lspec_2.skipcity = max([z modz],[],2) < rbspec.thresh;
        end
        lspec_2.skipcity(1) = 0;
    end
    
    % Standard adjustment
    adj = z2>rbspec.thresh &  lspec_2.skipcity(2:end)==0;
    smdevsz2(adj==1) =  smdevsz(adj==1) + 0.01*rbspec.step*w2(adj==1);
    
    if spec.modsm==1
        adjmod = modz>rbspec.thresh & lspec_2.skipcity==0;
        modsmz2(adjmod==1) =  modsmz(adjmod==1) - 0.01*rbspec.step*modw(adjmod==1);
    else
        modsmz2 = modsmz;
    end
    
    % Avoiding negative supermarkups
    smdevsz2(adj==1 & smdevsz2 < -smbase0) = 0.5*(smbase0+smdevsz2(adj==1 & smdevsz2 < -smbase0));
    smdevsz2(adj==2 & smdevsz2 < -smbase0) = - smbase0;
    
    % Resetting bad values: set mkt m to 5% of baseline
    smdevsz2(adj==1 & badval==1) = -0.95*smbase0;
    
    % Assessing balance for next step
    smvecz = [smbase0 ; smdevsz2];
    [~,~,~,balz2,~,~] = f_loss_bind(smvecz,mktvec,vars,ids,idat,daugfile,spec,lspec_2,[],modsmz);
           
    % Saving for next iteration
    smdevsz = smdevsz2;
    baly = balz2.rat;
    badval = balz2.num(2:end) < 0 ;
    badval2 = balz2.den(2:end) < 0 ;
    
    if spec.modsm==1
        modsmz = modsmz2;
        modbaly = balz2.ratmod;
    end
    
    
    if(max(badval(2:end))==1)
       disp('BADVAL PRESENT: dpi^PL / dm < 0'); 
       disp(badval(2:end)');
    end
    
    if(max(badval2(2:end))==1)
       disp('BADVAL PRESENT: dpi^PL > dpi^D'); 
       disp(badval2(2:end)');
    end
    
    % Progress Report
    if spec.modsm==0
        allrats = balz2.rat;
    elseif spec.modsm==1
        allrats = [balz2.rat  balz2.ratmod];
    end
    if spec.pool==1
        disp(allrats(lspec_2.skipcity==0,:)');
    elseif spec.pool==0
        disp(allrats');
    end
    
    % Convergence criterion
    if  max(max(abs(allrats-1))) < rbspec.thresh 
        tc = 9999;
    else
        tc = tc+1;
    end
    
end


% Resetting in order to get the final accounting (below)
if spec.pool==1 
    lspec_2.skipcity = zeros(size(lspec_2.skipcity));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the Rebalancing, pick the base sm to make ICC bind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Re-impute base value, using rebalanced deviations
lspec_2.fullvec = 0;
lspec_2.getder = 0;
fb = @(smvec)f_loss_bind(smvec,mktvec,vars,ids,idat,daugfile,spec,lspec_2,smdevsz,modsmz);
[smbase2,~] = fminsearch(fb,smbase0,spec.options.simp);

% Output objects
smout = [smbase2 ; smdevsz];
modsmout = modsmz;


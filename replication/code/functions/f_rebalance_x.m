function [smbase1,smdevs1,corner] = ...
    f_rebalance_x(smbase0,mktvec,vars,ids,idat,daugfile,spec,lspec,rbspec,smdevs0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the region-specific supermarkup deviations that
% satify the FOCs, given a base supermarkup. Use for case of size-specific
% supermarkups.
%
% Called by:
%   - main_supply_bind_x.m
% Calls the following user-specified functions:
%   - f_loss_bind_x.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% An initial assessment of balance
lspec_2 = lspec;
lspec_2.getder = 1;
xbase = smbase0*rbspec.scalebase;
[~,~,~,~,~,bal0,~] = f_loss_bind_x(xbase,mktvec,vars,ids,idat,daugfile,spec,lspec_2,smdevs0);

disp([ bal0.rat ]');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ad-Hoc Rebalancing Procedure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tic

balz = bal0.lev;
baly = bal0.rat;
badval = zeros(size(smdevs0));
badval2 = zeros(size(smdevs0));
smdevsz = smdevs0;
corner = 0;      
cornertest = [0 0 0]';

tc = 2;
while tc < 200    
  
    smdevsz2 = smdevsz;
    z = abs(baly-1); 
    %w = balz;
    w = baly-1;
    w(w>0.5)=0.5;
    w(w<-0.5)=-0.5;
    
    % Adjust city-skipper
    if spec.pool==1
        lspec_2.skipcity = max(z,[],2) < rbspec.thresh;
        lspec_2.skipcity(1) = 0;
    end
    
    % Adjust deviation
    adj = z>rbspec.thresh & lspec_2.skipcity==0;
    smdevsz2(adj==1) = smdevsz(adj==1) + 0.01*rbspec.step*w(adj==1);

    
    % Resetting bad values: set mkt m to 5% of baseline
    if size(baly,1)>1
        smdevsz2(adj==1 & badval==1) = -0.30*rbspec.scalebase*smbase0;
    else
        smdevsz2(adj==1 & baly<0) =  -0.30*rbspec.scalebase*smbase0;
    end
        
        
    % Assessing balance for next step
    [~,~,~,~,~,balz2,~] = f_loss_bind_x(xbase,mktvec,vars,ids,idat,daugfile,spec,lspec_2,smdevsz2);
    
    % Saving for next step
    smdevsz = smdevsz2;
    balz = balz2.lev;
    baly = balz2.rat;
    
    % Bad happenings
    if size(balz2.num,1)>1
        
        badval(2:end,:) = balz2.num(2:end,:) < 0 ;
        badval2(2:end,:) = balz2.den(2:end,:) < 0 ;
        
        % Alerts for bad happenings
        if(max(badval(2:end))==1)
            disp('BADVAL PRESENT: dpi^PL / dm < 0');
            %disp(badval(2:end)');
        end
        if(max(badval2(2:end))==1)
            disp('BADVAL PRESENT: dpi^PL > dpi^D');
            %disp(badval2(2:end)');
        end
        
    end
    
    % Progress Report
    temp = balz2.rat;
    if spec.pool==1
        disp(temp(lspec_2.skipcity==0,:)');
    elseif spec.pool==0
        disp(temp');
    end
        
    % Identifying corner solutions
    if spec.pool==0 && spec.bysize==2
       cornertest(3) = cornertest(2);
       cornertest(2) = cornertest(1);
       cornertest(1) = balz2.rat(2);
       if cornertest(3)>1 && cornertest(2)<1 && cornertest(1)>1
           corner = 1;
       elseif cornertest(3)<1 && cornertest(2)>1 && cornertest(1)<1
           corner = 1;
       end
    end
    
    % Convergence criterion
    tc = tc+1;
    if ( max(abs(balz2.rat-1)) < rbspec.thresh )
        tc = 9999;
    end
    if corner==1
        tc = 9999;
    end

end

%toc/60


% Balanced deviations (a function output argument)
smdevs1 = smdevsz;

% Resetting in order to get the final accounting (below)
if spec.pool==1 
    lspec_2.skipcity = zeros(size(smdevs1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the Rebalancing, pick the base sm to make ICC bind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tic

% Re-impute base value, using balanced deviations
lspec_2.getder = 0;
fb = @(smbase)f_loss_bind_x(smbase,mktvec,vars,ids,idat,daugfile,spec,lspec_2,smdevs1);
[smbase1,~] = fminsearch(fb,smbase0,lspec_2.optsimp);

%toc/60


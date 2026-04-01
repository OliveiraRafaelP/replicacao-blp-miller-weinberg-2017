function main_supply_bind(path,spec,fy)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function imputes supermarkups and marginal costs under assumption
% of **binding** ICCs.  The input arguments specify a discount factor 
% ("df"). The model features region-specific supermarkups with pooled ICCs.
%
% Called by:
%   - main_v7.m
% Calls the following user-specified functions:
%   - main_data.m
%   - f_loss_bind.m
%   - f_rebalance.m
% Saved results are loaded in:
%   - combine_imputed.m
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

%Pick city*periods in the fiscal year
mktvec = unique(ids.cdid(ids.fiscid==fy))';


% Starting values: sm for r=1 and R-1 deviaions
numdev = length(spec.city_in)-1;
if fy <2008
    startval = [1.0 + spec.df ; zeros(numdev,1)];
elseif fy>2008
    startval = [2.5 + spec.df ; zeros(numdev,1)];
    if spec.df >= 0.4
        startval = [3 + spec.df ; zeros(numdev,1)];
    end
end
clear numdev

% Initial imputation options
lspec.getder = 1;       % get the numerical derivatives
lspec.fullvec = 1;      % evaluate holding fixed deviations (0) or based on full sm vector (1) 
lspec.icc_restrict = 1; % restrict ICC analysis to Coors/MillerCoors
lspec.scale = 100;      % scale up balancing moments
lspec.purpose = 'I';    % imputation
lspec.mcleader = 0;     % ABI is leader, not MillerCoors

% This allows some cities to be skipped in rebalancing
lspec.skipcity = zeros(length(spec.city_in),1);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system of equations using the obtained starting value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration 0: hold fix deviations and make ICC bind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smdevs = startval(2:end);

smbase = startval(1);
lspec.fullvec = 0;
lspec.getder = 0;
fb = @(smvec)f_loss_bind(smvec,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs,[]);
[smbase0] = fminsearch(fb,smbase,spec.options.simp);

disp(smbase0)

if smbase0<0.05
    disp('Converged to m=0: Relaunch with New Starting Vale');
    BREAK
end


% Storing preliminary results

cd(strcat(path.data2,'/',spec.ssubfolder))

if fy==2006
    
    save sres_bind_2006 smbase0
    
elseif fy==2007
    
    save sres_bind_2007 smbase0
        
elseif fy==2010
    
    save sres_bind_2010 smbase0
    
elseif fy==2011
    
    save sres_bind_2011 smbase0
        
end

cd(strcat(path.code1))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebalance and then make the ICC bind again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smvec0 = [smbase0 ; smdevs];

% Iteration 1: 
rbspec.thresh = 0.2;
rbspec.step = 5.0;
[smout_1,bal_1] = f_rebalance(smvec0,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,[]);

disp('Iteration 1 Results:')
disp([bal_1.rat'])

% Iteration 2
rbspec.thresh = 0.1;
[smout_2,bal_2] = f_rebalance(smout_1,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,[]);

disp('Iteration 2 Results:')
disp([bal_2.rat'])


% Iteration 3 tightens criterion
rbspec.thresh = 0.05;
[smout_3,bal_3] = f_rebalance(smout_2,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,[]);

disp('Iteration 3 Results:')
disp([bal_3.rat'])

% Iteration 4 tightens criterion
rbspec.thresh = 0.025;
[smout_4,bal_4] = f_rebalance(smout_3,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,[]);

disp('Iteration 4 Results:')
disp([bal_4.rat'])


% Iteration 5 tightens criterion again
rbspec.thresh = 0.010;
[smout_5,bal_5] = f_rebalance(smout_4,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,[]);

disp('Iteration 5 Results:')
disp([bal_5.rat'])

% Iteration 6 tightens criterion again
rbspec.thresh = 0.005;
[smout_6,bal_6] = f_rebalance(smout_5,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,[]);

disp('Iteration 6 Results:')
disp([bal_6.rat'])

% Iteration 7 tightens criterion again
rbspec.thresh = 0.001;
[smout_7,bal_7] = f_rebalance(smout_6,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,[]);

disp('Iteration 7 Results:')
disp([bal_7.rat'])

% Iteration 8 repeats to confirm
[smout_8,bal_8] = f_rebalance(smout_7,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,[]);

disp('Iteration 8 Results:')
disp([bal_8.rat'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Economic staticstics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smfinal = smout_8;

%Imputing other statistics
lspec.fullvec = 1;
lspec.getder = 1;
lspec.icc_restrict = 0;
[~,~,outdata,balance,bindid,gicc] = f_loss_bind(smfinal,mktvec,vars,ids,1,daugfile,spec,lspec,[],[]);
        
mc = outdata.mc;
pnash = outdata.pnash;
snash = outdata.snash;

sdays = toc/60/60/24;

%disp(sdays);
%disp(df);
%disp(smfinal');

%%%%%%%%%%%%%%%%%%
% Storing results
%%%%%%%%%%%%%%%%%%

cd(strcat(path.data2,'/',spec.ssubfolder))

if fy==2006
    
    save sres_bind_2006 smbase0 mc pnash snash smfinal balance sdays bindid gicc outdata
    
elseif fy==2007
    
    save sres_bind_2007 smbase0 mc pnash snash smfinal balance sdays bindid gicc outdata
    
elseif fy==2010
    
    save sres_bind_2010 smbase0 mc pnash snash smfinal balance sdays bindid gicc outdata
    
elseif fy==2011
    
    save sres_bind_2011 smbase0 mc pnash snash smfinal balance sdays bindid gicc outdata
    
end

cd(strcat(path.code1))


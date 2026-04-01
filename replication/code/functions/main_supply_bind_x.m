function main_supply_bind_x(path,spec,fy)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function imputes supermarkups and marginal costs under assumption
% of **binding** ICCs.  The input arguments specify a discount factor 
% ("df"). The model features region-specific supermarkups with pooled ICCs.
%
% Called by:
%   - main_v6.m
% Calls the following user-specified functions:
%   - main_data.m
%   - f_loss_bind_x.m
%   - f_rebalance_x.m
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


spec.fy = fy;

% Pick city*periods in the fiscal year
mktvec = unique(ids.cdid(ids.fiscid==fy))';

% Starting value for the base supermarkup
if fy <2008
    smbase = 1.0 + spec.df;
elseif fy>2008
    smbase = 2.5 + spec.df;
    if spec.df >= 0.4
        smbase = 3 + spec.df;
    end
end
clear numdev

% Starting values for the deviations
if spec.bysize==0
    smdevs0 = zeros(length(spec.city_in),1);
elseif spec.bysize==1
    smdevs0 = zeros(length(spec.city_in),3);
elseif spec.bysize==2
    smdevs0 = zeros(length(spec.city_in),2);
end


% Optimization options: only use tight
lspec.optsimp = optimset('GradObj','off','MaxIter',50000,'MaxFunEvals',100000,...
    'Display','iter','TolFun',1e-4,'TolX',1e-4);

% Initial imputation options
lspec.getder = 0;       % get the numerical derivatives
lspec.icc_restrict = 1; % restrict ICC analysis to Coors/MillerCoors
lspec.nd2 = 1;          % 0=one-sided numerical derivatives (1=two-sided)
lspec.purpose = 'I';    % Imputing not simulating

% This allows some cities to be skipped in rebalancing
lspec.skipcity = zeros(size(smdevs0,1),1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system of equations using the obtained starting value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration 0: hold fix deviations and make ICC bind
%   - obtains optimal single supermarkup (across regions and sizes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fb = @(smbase)f_loss_bind_x(smbase,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs0);
[smbase0] = fminsearch(fb,smbase,lspec.optsimp);

disp(smbase0)

if smbase0<0.05
    disp('Converged to m=0: Relaunch with New Starting Vale');
    BREAK
end

% Storing preliminary results

cd(strcat(path.data2,'/',spec.ssubfolder))

if fy==2006
    
    save sres_bind_x_2006 smbase0
    
elseif fy==2007
    
    save sres_bind_x_2007 smbase0
    
elseif fy==2010
    
    save sres_bind_x_2010 smbase0
    
elseif fy==2011
    
    save sres_bind_x_2011 smbase0

end

cd(strcat(path.code1))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebalance and then make the ICC bind again
%   - Each iteration has a max of ~6 hours processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Iteration 1: 
rbspec.thresh = 0.1;
rbspec.step = 5.0;
rbspec.scalebase = 0.50;
lspec.nd2 = 1; 
[smbase1,smdevs1] = f_rebalance_x(smbase0,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,smdevs0);

disp('Iteration 1 Done')

% Iteration 1b: 
rbspec.thresh = 0.05;
rbspec.step = 5.0;
rbspec.scalebase = 0.70;
lspec.nd2 = 1; 
[smbase1b,smdevs1b] = f_rebalance_x(smbase1,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,smdevs1);

disp('Iteration 1b Done')

% Iteration 1c: 
rbspec.thresh = 0.05;
rbspec.step = 5.0;
rbspec.scalebase = 0.90;
lspec.nd2 = 1; 
[smbase1c,smdevs1c] = f_rebalance_x(smbase1b,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,smdevs1b);

disp('Iteration 1c Done')

% Iteration 2
rbspec.thresh = 0.05;
rbspec.step = 5.0;
rbspec.scalebase = 1.00;
[smbase2,smdevs2] = f_rebalance_x(smbase1c,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,smdevs1c);

disp('Iteration 2 Done')
 
% Iteration 3: 
rbspec.thresh = 0.01;
rbspec.step = 5;
[smbase3,smdevs3] = f_rebalance_x(smbase2,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,smdevs2);

disp('Iteration 3 Done')


% Iteration 4
rbspec.thresh = 0.001;
rbspec.step = 5;
[smbase4,smdevs4] = f_rebalance_x(smbase3,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,smdevs3);

disp('Iteration 4 Done')


% Iteration 5 repeats
rbspec.thresh = 0.001;
rbspec.step = 5;
[smbase5,smdevs5] = f_rebalance_x(smbase4,mktvec,vars,ids,1,daugfile,spec,lspec,rbspec,smdevs4);

disp('Iteration 5 Done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Economic statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smbase = smbase5;
smdevs = smdevs5;


lspec.icc_restrict = 0;
lspec.getder = 1;

[~,gbind,mc,pnash,snash,balance,bindid,gicc] = ...
    f_loss_bind_x(smbase,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs);

sdays = toc/60/60/24;


%%%%%%%%%%%%%%%%%%
% Storing results
%%%%%%%%%%%%%%%%%%

cd(strcat(path.data2,'/',spec.ssubfolder))

if fy==2006
    
    save sres_bind_x_2006 smbase0 mc pnash snash smbase smdevs balance sdays bindid gbind gicc
    
elseif fy==2007
    
    save sres_bind_x_2007 smbase0 mc pnash snash smbase smdevs balance sdays bindid gbind gicc
    
elseif fy==2010
    
    save sres_bind_x_2010 smbase0 mc pnash snash smbase smdevs balance sdays bindid gbind gicc
    
elseif fy==2011
    
    save sres_bind_x_2011 smbase0 mc pnash snash smbase smdevs balance sdays bindid gbind gicc
    
end

cd(strcat(path.code1))


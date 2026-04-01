function results_cfmergers(path,scen)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes equilibrium in the five merger simulations, as
% well as in the counterfactual featuring MillerCoors as a price leader.
%
% Called by:
%   - main_v7.m
% Calls the following user-specified functions:
%   - main_spec.m
%   - main_data.m
%   - f_loss_bind.m
%   - f_rebalance.m
%   - f_mu.m
%   - rcnl_inshr.m
%   - rcnl_der1.m
% Reads in data created in the following functions:
%   - combine_imputed.m
%   - results_costregs.m
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
pcoefi = daugfile.pcoefi(ids.obsindemand==1,:);
deltanp = daugfile.deltanp(ids.obsindemand==1);
rho = daugfile.rho;

% Results at selected discount factor
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc pnash snash sm_07 sm_10 sm_11 ;
cd(path.code1)

% Loading cost parameters
cd(strcat(path.data2))
load gammapar gamma 
cd(path.code1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing data for simulations, scenario-specific adjustments to
% ownership, marginal costs, time periods, and starting values for m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if scen==1  % Unwind Miller/Coors (ownership and efficiencies)
    
    % Will evaluate FY 2010
    mktvec = unique(ids.cdid(ids.fiscid==2010))';
    
    % Changing ownership/coalition structure
    idx = ids;
    idx.firmid(idx.brndid==3) = 4;
    idx.firmid(idx.brndid==4) = 4;
    idx.coalid = (idx.firmid==1|idx.firmid==5|idx.firmid==4);
    
    % Changing marginal costs
    mcx = mc;
    mcx(idx.firmid==5) = mcx(idx.firmid==5) - gamma(2);    
    mcx(idx.firmid==4) = mcx(idx.firmid==4) - gamma(3);
    mcx = mcx - gamma(4)*vars.dist + gamma(4)*vars.distbutfor;
    
    % Starting value for supermarkups
    smdevs0 = sm_07(2:end) - sm_07(1);
    smbase00 = sm_07(1);
    
elseif scen==2  % Remove Miller/Coors efficiencies, keep merger
    
    % Will evaluate FY 2010
    mktvec = unique(ids.cdid(ids.fiscid==2010))';
    
    % Unchanged ownership/coalition
    idx = ids;
    
    % Changing marginal costs
    mcx = mc;
    mbrands = idx.brndid==11 | idx.brndid==12 | idx.brndid==13;
    cbrands = idx.brndid== 3 | idx.brndid== 4;
    mcx(mbrands==1) = mcx(mbrands==1) - gamma(2);
    mcx(cbrands==1) = mcx(cbrands==1) - gamma(3);
    mcx = mcx - gamma(4)*vars.dist + gamma(4)*vars.distbutfor;
    clear mbrands cbrands
    
    % Starting value for supermarkups
    smdevs0 = 0.0*(sm_10(2:end) - sm_10(1));
    smbase00 = sm_10(1);
    
    
elseif scen==6  % ABI/Modelo no efficiencies
    
    % Will evaluate FY 2011
    mktvec = unique(ids.cdid(ids.fiscid==2011))';
    
    % Changing ownership/coalition structure
    idx=ids;
    idx.firmid(idx.brndid==5) = 1;
    idx.firmid(idx.brndid==6) = 1;
    idx.coalid = (idx.firmid==1|idx.firmid==5|idx.firmid==4);
    
    % Unchanged marginal costs
    mcx = mc;
    
    % Starting value for supermarkups
    smdevs0 = sm_11(2:end) - sm_11(1);
    smbase00 = 0.8*sm_11(1);
    
elseif scen==7  % ABI/Modelo minor efficiencies
   
    % Will evaluate FY 2011
    mktvec = unique(ids.cdid(ids.fiscid==2011))';
    
    % Changing ownership/coalition structure
    idx=ids;
    idx.firmid(idx.brndid==5) = 1;
    idx.firmid(idx.brndid==6) = 1;
    idx.coalid = (idx.firmid==1|idx.firmid==5|idx.firmid==4);
    
    % Minor efficiencies
    mcx = mc;
    mcx = mcx - 0.50*(ids.firmid==2);
    
    % Starting value for supermarkups
    smdevs0 = sm_11(2:end) - sm_11(1);
    smbase00 = 0.8*sm_11(1);
    
elseif scen==8 % ABI/Modelo major efficiencies
    
    % Will evaluate FY 2011
    mktvec = unique(ids.cdid(ids.fiscid==2011))';
    
    % Changing ownership/coalition structure
    idx = ids;
    idx.firmid(idx.brndid==5) = 1;
    idx.firmid(idx.brndid==6) = 1;
    idx.coalid = (idx.firmid==1|idx.firmid==5|idx.firmid==4);
    
    % Major efficiencies---no Bertrand price increase---a number of steps
    mcx1 = mc;
    mcx2 = mc;
    
    % Adjusted mean valuations --- evaluate at Nash prices
    delta2 = deltanp + daugfile.alpha*pnash;

    % Looping through markets
    for mktid=mktvec
        
        owner1   = f_ownMat(ids.firmid(idx.cdid==mktid));
        owner2   = f_ownMat(idx.firmid(idx.cdid==mktid));

        pcoefiM  = pcoefi(ids.cdid==mktid,:);
        delta2M  = delta2(idx.cdid==mktid);
        x2M      = vars.x2(idx.cdid==mktid,:);
        x2M(:,1) = pnash(idx.cdid==mktid);
        dfullM   = vars.dfull(idx.cdid==mktid,:);
        
        [mu2M,~] = f_mu(daugfile.theta2w,vars.ns,dfullM,x2M);
        
        [shareiM,scondiM,sgroupiM]=rcnl_indsh(exp(delta2M),exp(mu2M),rho,0,0,1);
        der = rcnl_der1(pcoefiM,shareiM,scondiM,sgroupiM,daugfile.rho);
        
        markup1 = -(owner1.*(der'))\snash(ids.cdid==mktid);
        markup2 = -(owner2.*(der'))\snash(ids.cdid==mktid);
        mcx1(ids.cdid==mktid) = pnash(ids.cdid==mktid) - markup1;
        mcx2(ids.cdid==mktid) = pnash(ids.cdid==mktid) - markup2;
        
    end
    
    mceff = mcx2-mcx1;
    mcx = mc + mceff;
    
    majoreff = [mean(mceff(ids.fiscid==2011 & ids.firmid==1)) ...
        mean(mceff(ids.fiscid==2011 & ids.firmid==2)) ];
        
    %Writing to .txt file---these go into footnote
    cd(strcat(path.figr))
        TT = table(round(majoreff,2));
        writetable(TT,'majoreff.txt','Delimiter',';');
    cd(path.code1)
    
    
    % Starting value for supermarkups
    smdevs0 = sm_11(2:end) - sm_11(1);
    smbase00 = 0.8*sm_11(1);
    
elseif scen==12 % MillerCoors is the leader---2010
    
    % Will evaluate FY 2010
    mktvec = unique(ids.cdid(ids.fiscid==2010))';
    
    % Same ownership/coalition
    idx=ids;
    
    % Same costs
    mcx = mc;
    
    % Starting values for supermarkups
    smdevs0 = 1.0*(sm_11(2:end) - sm_11(1));
    smbase00 = sm_11(1);
    
end

% Starting value for Modelo supermarkup
if scen==6 || scen==7 || scen==8
    spec.modsm = 1;
    modsm0 = zeros(size(spec.city_in))';
else
    modsm0 = [];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating Bertrand prices and shares, which are inputs to the PLE
% simulation code, as it makes sense to do this calculation only once.
%   - Takes a couple minutes to solve each market.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pnashx = pnash;
snashx = snash;

% Looping through markets
for mktid=mktvec
    
    x2M = vars.x2(idx.cdid==mktid,:);
    pcoefiM = pcoefi(idx.cdid==mktid,:);
    deltanpM = deltanp(idx.cdid==mktid);
    dfullM = vars.dfull(idx.cdid==mktid,:);    
    mcxM = mcx(idx.cdid==mktid);
    owner   = f_ownMat(idx.firmid(idx.cdid==mktid));

    p = pnash(idx.cdid==mktid);
    popt = ones(size(p));
    pfix = [];
    
    f = @(p)cf_foc_partial(p,daugfile.alpha,daugfile.rho,deltanpM,mcxM,...
        owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns);
    [fpn] = fsolve(f,p,spec.options.cf);
    
    [~,fsn,~,~] = cf_foc_partial(fpn,daugfile.alpha,daugfile.rho,deltanpM,mcxM,...
        owner,pfix,popt,x2M,pcoefiM,dfullM,daugfile.theta2w,vars.ns); 
    
    pnashx(idx.cdid==mktid) = fpn;
    snashx(idx.cdid==mktid) = fsn;
    
    clear owner p popt pfix fpn
         
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conducting the PLE simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting organized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% Input data
idat.mc = mcx;
idat.pnash = pnashx;
idat.snash = snashx;

% Simulation options
lspec.getder = 0;       
lspec.fullvec = 0;      
lspec.icc_restrict = 0; % restrict ICC analysis to Coors/MillerCoors
lspec.scale = 100;      
lspec.purpose = 'S';    
lspec.skipcity = zeros(length(spec.city_in),1);
lspec.mcleader = 0;

% This modification is needed to use the skipcity functionality
% if scen<=8
%     lspec.icc_restrict = 0; % restrict ICC analysis to Coors/MillerCoors for computational savings
% end

% This is the key line for spec 12
if scen==12
    lspec.mcleader = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Adjusting supermarkup levels (holding fixed deviations)
fb = @(smvec)f_loss_bind(smvec,mktvec,vars,idx,idat,daugfile,spec,lspec,smdevs0,modsm0);
[smbase0,~,~] = fminsearch(fb,smbase00,spec.options.simp);

if smbase0<0.05
    disp('Converged to m=0: Relaunch with New Starting Value');
    BREAK
end

smvec0 = [smbase0 ; smdevs0];


% Iteration 1
rbspec.thresh = 0.1;
rbspec.step = 5.0;
[smout_1,bal_1,modsm_1] = f_rebalance(smvec0,mktvec,vars,idx,idat,daugfile,spec,lspec,rbspec,modsm0);

%f_rebalance(smvec0,mktvec,vars,ids,idat,daugfile,spec,lspec,rbspec,modsm0)


disp('Iteration 1 Results:')
disp([bal_1.rat'])

% Different treatment for scenario 7 with the fake data due to numerics
if scen==7
    
    % Iteration 1b
    rbspec.thresh = 0.05;
    rbspec.step = 1.0;
    [smout_1b,bal_1b,modsm_1b] = f_rebalance(smout_1,mktvec,vars,idx,idat,daugfile,spec,lspec,rbspec,modsm_1);
    
    smfinal = smout_1b;
    modsmfinal = modsm_1b;
    
else
        
    % Iteration 2
    rbspec.thresh = 0.01;
    [smout_2,bal_2,modsm_2] = f_rebalance(smout_1,mktvec,vars,idx,idat,daugfile,spec,lspec,rbspec,modsm_1);
    
    disp('Iteration 2 Results:')
    disp([bal_2.rat'])
    
    
    
    %Iteration 3 --- smaller steps
    rbspec.thresh = 0.001;
    rbspec.step = 2.5;
    [smout_3,bal_3,modsm_3] = f_rebalance(smout_2,mktvec,vars,idx,idat,daugfile,spec,lspec,rbspec,modsm_2);
    
    disp('Iteration 3 Results:')
    disp([bal_3.rat'])
    
    % Iteration 4
    rbspec.thresh = 0.001;
    [smout_4,bal_4,modsm_4] = f_rebalance(smout_3,mktvec,vars,idx,idat,daugfile,spec,lspec,rbspec,modsm_3);
    
    disp('Iteration 4 Results:')
    disp([bal_4.rat'])
    
    % Iteration 5
    rbspec.thresh = 0.001;
    [smout_5,bal_5,modsm_5] = f_rebalance(smout_4,mktvec,vars,idx,idat,daugfile,spec,lspec,rbspec,modsm_4);
    
    disp('Iteration 5 Results:')
    disp([bal_5.rat'])
    
    smfinal = smout_5;
    modsmfinal = modsm_5;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Economic statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



sm_cf = smfinal;
sm_cf(2:end) = sm_cf(1)+sm_cf(2:end);

%Imputing other statistics
lspec.fullvec = 1;
lspec.getder = 1;
lspec.icc_restrict = 0;
[~,~,outdata,balance,bindid,gicc] = f_loss_bind(smfinal,mktvec,vars,idx,idat,daugfile,spec,lspec,[],modsmfinal);

sdays = toc/60/60/24;


%%%%%%%%%%%%%%%%%%
% Storing results
%%%%%%%%%%%%%%%%%%

cd(strcat(path.data2))

if scen==1
    
    save cfscen_01 outdata smfinal balance sdays bindid gicc sm_cf
    
elseif scen==2
    
    save cfscen_02 outdata smfinal balance sdays bindid gicc sm_cf
        
elseif scen==6
    
    save cfscen_06 outdata smfinal balance sdays bindid gicc sm_cf modsmfinal
    
elseif scen==7
    
    save cfscen_07 outdata smfinal balance sdays bindid gicc sm_cf modsmfinal
    
elseif scen==8
    
    save cfscen_08 outdata smfinal balance sdays bindid gicc sm_cf modsmfinal
    
elseif scen==12
    
    save cfscen_12 outdata smfinal balance sdays bindid gicc sm_cf
    
end


cd(strcat(path.code1))

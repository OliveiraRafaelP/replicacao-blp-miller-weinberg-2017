function results_baseanalysis(path)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces tables and figures that appear in the text, along
% with associated numbers.
%
% Called by:
%   - main_v7.m
% Calls the following user-specified functions:
%   - main_spec.m
%   - main_data.m
%   - f_loss_bind.m
%   - f_inclusive.m
% Reads in data created in the following functions:
%   - combine_imputed.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Baseline specification
dfx = 26;
df = dfx/100;
spec = main_spec(df);
spec.bysize = 0;  

% Grabbing data and defining fiscnum
[vars,ids] = main_data(path,spec);
fiscnum = grp2idx(ids.fiscid);

% Demand results
cd(strcat(path.data2))
load daugfile daugfile;
cd(path.code1)
pcoefi = daugfile.pcoefi;
deltanp = daugfile.deltanp;
mu = daugfile.mu;
ai = daugfile.ai;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 2: average price and quantities (conditional volume shares)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condshr = vars.s_jt ./ vars.inshr;

% Table construction
out = zeros(length(unique(ids.brndid)),7);
n=1;
for i = unique(ids.brndid)'
    s_06 = mean(condshr(ids.fiscid==2011 & ids.brndid==i & ids.sizeid==1));
    p_06 = mean(vars.p_jt(ids.fiscid==2011 & ids.brndid==i & ids.sizeid==1));
    s_12 = mean(condshr(ids.fiscid==2011 & ids.brndid==i & ids.sizeid==2));
    p_12 = mean(vars.p_jt(ids.fiscid==2011 & ids.brndid==i & ids.sizeid==2));
    s_24 = mean(condshr(ids.fiscid==2011 & ids.brndid==i & ids.sizeid==3));
    p_24 = mean(vars.p_jt(ids.fiscid==2011 & ids.brndid==i & ids.sizeid==3));
    s_tot = s_06 + s_12 + s_24;
    out(n,:) = [round(s_06,3) round(p_06,2) ...
                round(s_12,3) round(p_12,2) ...
                round(s_24,3) round(p_24,2) ...
                round(s_tot,3)];
    n=n+1;
end


%Writing to .txt file
cd(strcat(path.figr))
    TT = table(out);
    writetable(TT,'sumstats_pq.txt','Delimiter',';');
cd(path.code1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 4.  Summary statistics by year:
%   - Mean supermarkup (average across regions)
%   - % Change in Total Profit
%   - Change in CS / Change in Profit
%   - % Change in Profit--by firm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Results at selected discount factor
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc pnash snash;
cd(path.code1)

% Shell file 
resmat = zeros(9,4);

% Overall Profit Data
prof_ple  = (vars.p_jt-mc).*vars.s_jt.*vars.msize;
prof_nash = (pnash-mc).*snash.*vars.msize;
resmat(2,:) = (accumarray(fiscnum,prof_ple) ./ accumarray(fiscnum,prof_nash)) - 1;


% Profit by firm
for f = unique(ids.firmid)'
    temp1 = prof_ple .* (ids.firmid==f);
    temp2 = prof_nash .* (ids.firmid==f);
    res = (accumarray(fiscnum,temp1) ./ accumarray(fiscnum,temp2)) - 1;
    if f==1
        resmat(3,:) = res;
    elseif f==5
        resmat(4,3:4) = res(3:4);
        resmat(5,1:2) = res(1:2);
    elseif f==4
        resmat(6,1:2) = res(1:2);
    elseif f==2
        resmat(7,:) = res;
    elseif f==3
        resmat(8,:) = res;
    end
    
end

% Changing profit units
resmat = resmat*100;

deltanpx = deltanp(ids.obsindemand==1);
mux = mu(ids.obsindemand==1,:);
aix = ai(ids.obsindemand==1,:);

% Consumer Surplus
delta_ple = deltanpx+daugfile.alpha*vars.p_jt;
delta_nash = deltanpx+daugfile.alpha*pnash;
mu_ple = mux;
mu_nash = mux-aix.*repmat(vars.p_jt,[1 size(aix,2)])+aix.*repmat(pnash,[1 size(aix,2)]);
cs_ple = f_inclusive(delta_ple,mu_ple,ids.cdindex,pcoefi,daugfile.rho,vars.msize);
cs_nash = f_inclusive(delta_nash,mu_nash,ids.cdindex,pcoefi,daugfile.rho,vars.msize);
dcs = accumarray(fiscnum(ids.cdindex),cs_nash) - accumarray(fiscnum(ids.cdindex),cs_ple);
resmat(end,:) = - dcs' ./ (accumarray(fiscnum,prof_ple) - (accumarray(fiscnum,prof_nash)))';

resmat = resmat([2;end],:);

%Writing to .txt file
cd(strcat(path.figr))
    TT1 = table(round(resmat,2));
    writetable(TT1,'eqeffects.txt','Delimiter',';');
cd(path.code1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 5: brewer markups
%   - also assorted markup statistics quoted in text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Results at selected discount factor
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc;
cd(path.code1)

% Markups
markups = vars.p_jt - mc;

stat1 = mean(markups);
stat2 = mean(vars.p_jt);
stat3 = mean(markups./vars.p_jt);
stat4 = stat2*0.65;
stat5 = mean(markups./stat4);
stat6 = stat2*0.35;

% Table construction
out = zeros(length(unique(ids.brndid)),6);
n=1;
for i = unique(ids.brndid)'
    pre06 = mean(markups(ids.fiscid==2007 & ids.brndid==i & ids.sizeid==1));
    pos06 = mean(markups(ids.fiscid==2010 & ids.brndid==i & ids.sizeid==1));
    pre12 = mean(markups(ids.fiscid==2007 & ids.brndid==i & ids.sizeid==2));
    pos12 = mean(markups(ids.fiscid==2010 & ids.brndid==i & ids.sizeid==2));
    pre24 = mean(markups(ids.fiscid==2007 & ids.brndid==i & ids.sizeid==3));
    pos24 = mean(markups(ids.fiscid==2010 & ids.brndid==i & ids.sizeid==3));    
    out(n,:) = [pre06 pos06 pre12 pos12 pre24 pos24];
    n=n+1;
end


%Writing to .txt file
cd(strcat(path.figr))
    TT = table(round(out,2));
    TT2 = table([stat1 stat2 stat3 stat4 stat5 stat6]);
    writetable(TT,'mean_markups.txt','Delimiter',';');
    writetable(TT2,'markup_stats.txt','Delimiter',';');
cd(path.code1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4: Line graphs that trace out the slack functions
%   - 10 min/call * 4 fy * 150 steps * (1/60) mins/hr = 100 hrs
%   - Will do each FY separately to save intermediate results
%   - Consider only 2007 and 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load results at selected discount factor
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc sm_06 sm_07 sm_10 sm_11 pnash snash;
cd(path.code1)

% Specification for loss function
lspec.skipcity = zeros(37,1);
lspec.fullvec = 2;
lspec.getder = 0;
lspec.scale = 100;
lspec.purpose = 'S';
lspec.icc_restrict = 0;
lspec.mcleader = 0;

% Other materials for loss function
smdevs = [];
idat.mc=mc;
idat.pnash=pnash;
idat.snash=snash;

% Number of steps
steps=150;

% Shell files
giccmat_06 = zeros(3,steps+1);
giccmat_07 = zeros(3,steps+1);
giccmat_10 = zeros(2,steps+1);
giccmat_11 = zeros(2,steps+1);
estats_06  = zeros(length(ids.cdid(ids.fiscid==2006)),9,steps+1); 
estats_07  = zeros(length(ids.cdid(ids.fiscid==2007)),9,steps+1); 
estats_10  = zeros(length(ids.cdid(ids.fiscid==2010)),9,steps+1); 
estats_11  = zeros(length(ids.cdid(ids.fiscid==2011)),9,steps+1); 


% FY 2006!
disp('Starting 2006');
mktvec = unique(ids.cdid(ids.fiscid==2006))';
for i=0:steps
    disp(i);
    [~,~,od,~,~,gicc] = f_loss_bind(sm_06*(i/100),mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,[]);
    giccmat_06(:,i+1) = gicc;
    estats_06(:,:,i+1) = [od.pnash od.snash od.pinash od.pcoll od.scoll od.picoll od.pdev od.sdev od.pidev];
end
cd(strcat(path.data2))
save baseslack_06 giccmat_06 estats_06
cd(strcat(path.code1))

% FY 2007!
disp('Starting 2007');
mktvec = unique(ids.cdid(ids.fiscid==2007))';
for i=0:steps
    disp(i);
    [~,~,od,~,~,gicc] = f_loss_bind(sm_07*(i/100),mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,[]);
    giccmat_07(:,i+1) = gicc;
    estats_07(:,:,i+1) = [od.pnash od.snash od.pinash od.pcoll od.scoll od.picoll od.pdev od.sdev od.pidev];
end
cd(strcat(path.data2))
save baseslack_07 giccmat_07 estats_07 
cd(strcat(path.code1))

% FY 2010!
disp('Starting 2010');
mktvec = unique(ids.cdid(ids.fiscid==2010))';
for i=0:steps
    disp(i);
    [~,~,od,~,~,gicc] = f_loss_bind(sm_10*(i/100),mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,[]);
    giccmat_10(:,i+1) = gicc;
    estats_10(:,:,i+1) = [od.pnash od.snash od.pinash od.pcoll od.scoll od.picoll od.pdev od.sdev od.pidev];
end
cd(strcat(path.data2))
save baseslack_10 giccmat_10 estats_10
cd(strcat(path.code1))

%FY 2011!
disp('Starting 2011');
mktvec = unique(ids.cdid(ids.fiscid==2011))';
for i=0:steps
    disp(i);
    [~,~,od,~,~,gicc] = f_loss_bind(sm_11*(i/100),mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,[]);
    giccmat_11(:,i+1) = gicc;
    estats_11(:,:,i+1) = [od.pnash od.snash od.pinash od.pcoll od.scoll od.picoll od.pdev od.sdev od.pidev];
end
cd(strcat(path.data2))
save baseslack_11 giccmat_11 estats_11
cd(strcat(path.code1))

% Loading data for graphing
cd(strcat(path.data2))
load baseslack_06 giccmat_06 
load baseslack_07 giccmat_07 
load baseslack_10 giccmat_10 
load baseslack_11 giccmat_11 
cd(strcat(path.code1))

% Prepping
smind = (0:150)/100;
bottom_06 = max(min(min(giccmat_06)),-max(max(giccmat_06)));
bottom_07 = max(min(min(giccmat_07)),-max(max(giccmat_07)));
bottom_10 = max(min(min(giccmat_10)),-max(max(giccmat_10)));
bottom_11 = max(min(min(giccmat_11)),-max(max(giccmat_11)));



% 2006 Slack Functions
figure
plot(smind*mean(sm_06),giccmat_06(1,:)', 'k', ...
    smind*mean(sm_06),giccmat_06(3,:)', 'k-.', ...
    smind*mean(sm_06),giccmat_06(2,:), 'k--')
set(gca,'YLim',[bottom_06,1.10*max(max(giccmat_07))])
set(gca,'XLim',[0,1.50*mean(sm_11)])
refline(1,0);
line([mean(sm_06) mean(sm_06)], get(gca, 'ylim'));
grid on
legend('ABI','Miller','Coors', 'Location', 'southwest','FontSize',18)
ylabel('Slack in IC Constraint','FontSize',20)
%xlabel('Supermarkup (Index)','FontSize',20)
xlabel('')
title('2006','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('iccobs_06'), 'pdf')
cd(path.code1)

% 2007 Slack Functions
figure
plot(smind*mean(sm_07),giccmat_07(1,:)', 'k', ...
    smind*mean(sm_07),giccmat_07(3,:)', 'k-.', ...
    smind*mean(sm_07),giccmat_07(2,:), 'k--')
set(gca,'XLim',[0,1.50*mean(sm_11)])
set(gca,'YLim',[bottom_07,1.10*max(max(giccmat_07))])
refline(1,0);
line([mean(sm_07) mean(sm_07)], get(gca, 'ylim'));
grid on
legend('ABI','Miller','Coors', 'Location', 'southwest','FontSize',18)
ylabel('Slack in IC Constraint','FontSize',20)
%ylabel('')
xlabel('Supermarkup','FontSize',20)
%xlabel('')
title('2007','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('iccobs_07'), 'pdf')
cd(path.code1)

% 2010 Slack Functions
figure
plot(smind*mean(sm_10),giccmat_10(1,:)', 'k', ...
    smind*mean(sm_10),giccmat_10(2,:), 'k--')
set(gca,'XLim',[0,1.50*mean(sm_11)])
set(gca,'YLim',[bottom_10,1.10*max(max(giccmat_10))])
refline(1,0);
line([mean(sm_10) mean(sm_10)], get(gca, 'ylim'));
grid on
legend('ABI','MillerCoors', 'Location', 'southwest','FontSize',18)
xlabel('Supermarkup','FontSize',20)
title('2010','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('iccobs_10'), 'pdf')
cd(path.code1)

% 2011 Slack Functions
figure
plot(smind*mean(sm_11),giccmat_11(1,:)', 'k', ...
    smind*mean(sm_11),giccmat_11(2,:), 'k--')
set(gca,'XLim',[0,1.50*mean(sm_11)])
set(gca,'YLim',[bottom_11,1.10*max(max(giccmat_10))])
refline(1,0);
line([mean(sm_11) mean(sm_11)], get(gca, 'ylim'));
grid on
legend('ABI','MillerCoors', 'Location', 'southwest','FontSize',18)
%ylabel('Slack in IC Constraint','FontSize',20)
ylabel('')
xlabel('Supermarkup','FontSize',20)
title('2011','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('iccobs_11'), 'pdf')
cd(path.code1)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3: Line graphs that trace out the profit functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading the data
cd(strcat(path.data2))
load baseslack_06 estats_06
load baseslack_07 estats_07
load baseslack_10 estats_10
load baseslack_11 estats_11
cd(strcat(path.code1))

smind = (0:150)/100;

% This provides the ordering of objects within estats_YY
%estats_07(:,:,i+1) = [od.pnash od.snash od.pinash od.pcoll od.scoll od.picoll od.pdev od.sdev od.pidev];

%%%% 2007 Profit terms %%%%

pi_nash_abi_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==1,3,:),[1 3 2]));
pi_nash_abi_07 = pi_nash_abi_07(end,:);

pi_nash_mill_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==5,3,:),[1 3 2]));
pi_nash_mill_07 = pi_nash_mill_07(end,:);

pi_nash_coor_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==4,3,:),[1 3 2]));
pi_nash_coor_07 = pi_nash_coor_07(end,:);

pi_coll_abi_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==1,6,:),[1 3 2]));
pi_coll_abi_07 = pi_coll_abi_07(end,:) ./ pi_nash_abi_07;

pi_coll_mill_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==5,6,:),[1 3 2]));
pi_coll_mill_07 = pi_coll_mill_07(end,:)  ./ pi_nash_mill_07;

pi_coll_coor_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==4,6,:),[1 3 2]));
pi_coll_coor_07 = pi_coll_coor_07(end,:) ./ pi_nash_coor_07;

pi_dev_abi_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==1,9,:),[1 3 2]));
pi_dev_abi_07 = pi_dev_abi_07(end,:) ./ pi_nash_abi_07;

pi_dev_mill_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==5,9,:),[1 3 2]));
pi_dev_mill_07 = pi_dev_mill_07(end,:) ./ pi_nash_mill_07;

pi_dev_coor_07 = cumsum(permute(estats_07(ids.firmid(ids.fiscid==2007)==4,9,:),[1 3 2]));
pi_dev_coor_07 = pi_dev_coor_07(end,:) ./ pi_nash_coor_07;


%%%% 2010 Profit terms %%%%

pi_nash_abi_10 = cumsum(permute(estats_10(ids.firmid(ids.fiscid==2010)==1,3,:),[1 3 2]));
pi_nash_abi_10 = pi_nash_abi_10(end,:);

pi_nash_mc_10 = cumsum(permute(estats_10(ids.firmid(ids.fiscid==2010)==5,3,:),[1 3 2]));
pi_nash_mc_10 = pi_nash_mc_10(end,:);

pi_coll_abi_10 = cumsum(permute(estats_10(ids.firmid(ids.fiscid==2010)==1,6,:),[1 3 2]));
pi_coll_abi_10 = pi_coll_abi_10(end,:) ./ pi_nash_abi_10;

pi_coll_mc_10 = cumsum(permute(estats_10(ids.firmid(ids.fiscid==2010)==5,6,:),[1 3 2]));
pi_coll_mc_10 = pi_coll_mc_10(end,:)  ./ pi_nash_mc_10;

pi_dev_abi_10 = cumsum(permute(estats_10(ids.firmid(ids.fiscid==2010)==1,9,:),[1 3 2]));
pi_dev_abi_10 = pi_dev_abi_10(end,:) ./ pi_nash_abi_10;

pi_dev_mc_10 = cumsum(permute(estats_10(ids.firmid(ids.fiscid==2010)==5,9,:),[1 3 2]));
pi_dev_mc_10 = pi_dev_mc_10(end,:) ./ pi_nash_mc_10;


% Profit figure: ABI 2007
figure
plot(smind*mean(sm_07),pi_coll_abi_07', 'k', ...
    smind*mean(sm_07),pi_dev_abi_07, 'k--')
    set(gca,'YLim',[1,1.50])
    line([mean(sm_07) mean(sm_07)], get(gca, 'ylim'));
    grid on
    legend('Price Leadership','Deviation', 'Location', 'southeast','FontSize',18)
    ylabel('Index Relative to Bertrand','FontSize',20)
    xlabel('Supermarkup','FontSize',20)
    title('ABI in 2007','FontSize',20)
    set(gcf,'color','white');
    set(gcf, 'PaperPosition', [0 0 5 5]); 
    set(gcf, 'PaperSize', [5 5]); 
cd(path.figr)
saveas(gcf, strcat('prof_abi_07'), 'pdf')
cd(path.code1)

% Profit figure: ABI 2010
figure
plot(smind*mean(sm_10),pi_coll_abi_10', 'k', ...
    smind*mean(sm_10),pi_dev_abi_10, 'k--')
    set(gca,'YLim',[1,1.50])
    line([mean(sm_10) mean(sm_10)], get(gca, 'ylim'));
    grid on
    legend('Price Leadership','Deviation', 'Location', 'southeast','FontSize',18)
    ylabel('Index Relative to Bertrand','FontSize',20)
    xlabel('Supermarkup','FontSize',20)
    title('ABI in 2010','FontSize',20)
    set(gcf,'color','white');
    set(gcf, 'PaperPosition', [0 0 5 5]); 
    set(gcf, 'PaperSize', [5 5]); 
cd(path.figr)
saveas(gcf, strcat('prof_abi_10'), 'pdf')
cd(path.code1)



% Profit figure: Miller 2007
figure
plot(smind*mean(sm_07),pi_coll_mill_07', 'k', ...
    smind*mean(sm_07),pi_dev_mill_07, 'k--')
    set(gca,'YLim',[1,1.50])
    line([mean(sm_07) mean(sm_07)], get(gca, 'ylim'));
    grid on
    legend('Price Leadership','Deviation', 'Location', 'southeast','FontSize',18)
    ylabel('')
    xlabel('Supermarkup','FontSize',20)
    title('Miller in 2007','FontSize',20)
    set(gcf,'color','white');
    set(gcf, 'PaperPosition', [0 0 5 5]); 
    set(gcf, 'PaperSize', [5 5]); 
cd(path.figr)
saveas(gcf, strcat('prof_mill_07'), 'pdf')
cd(path.code1)

% Profit figure: Coors 2007
figure
plot(smind*mean(sm_07),pi_coll_coor_07', 'k', ...
    smind*mean(sm_07),pi_dev_coor_07, 'k--')
    set(gca,'YLim',[1,1.50])
    line([mean(sm_07) mean(sm_07)], get(gca, 'ylim'));
    grid on
    legend('Price Leadership','Deviation', 'Location', 'southeast','FontSize',18)
    %ylabel('Index Relative to Bertrand','FontSize',20)
    ylabel('')
    xlabel('Supermarkup','FontSize',20)
    title('Coors in 2007','FontSize',20)
    set(gcf,'color','white');
    set(gcf, 'PaperPosition', [0 0 5 5]); 
    set(gcf, 'PaperSize', [5 5]); 
cd(path.figr)
saveas(gcf, strcat('prof_coor_07'), 'pdf')
cd(path.code1)

% Profit figure: MillerCoors 2010
figure
plot(smind*mean(sm_10),pi_coll_mc_10', 'k', ...
     smind*mean(sm_10),pi_dev_mc_10, 'k--')
    set(gca,'YLim',[1,1.50])
    line([mean(sm_10) mean(sm_10)], get(gca, 'ylim'));
    grid on
    legend('Price Leadership','Deviation', 'Location', 'southeast','FontSize',18)
    %ylabel('Index Relative to Bertrand','FontSize',20)
    ylabel('')
    xlabel('Supermarkup','FontSize',20)
    title('MillerCoors in 2010','FontSize',20)
    set(gcf,'color','white');
    set(gcf, 'PaperPosition', [0 0 5 5]); 
    set(gcf, 'PaperSize', [5 5]); 
cd(path.figr)
saveas(gcf, strcat('prof_mc_10'), 'pdf')
cd(path.code1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5: Scatterplots of m on shares and HHI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Getting variables: price, share, market size
spec = main_spec(0.26);
spec.bysize = 0;
[vars,ids] = main_data(path,spec);

% Load results at selected discount factor
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine sm_07 sm_10;
cd(path.code1)


% Firm shares by region
sharemat_07 = zeros(length(unique(ids.cityid)),5);
sharemat_10 = zeros(length(unique(ids.cityid)),5);
for f = 1:5
    selector_07 = (ids.fiscid==2007) .* (ids.firmid==f);
    selector_10 = (ids.fiscid==2010) .* (ids.firmid==f);
    sharemat_07(:,f)=accumarray(ids.cityid(selector_07==1),vars.s_jt(selector_07==1));
    if f ~= 4
        sharemat_10(:,f)=accumarray(ids.cityid(selector_10==1),vars.s_jt(selector_10==1));
    end
end
sharemat_07 = sharemat_07(spec.city_in,:);
sharemat_10 = sharemat_10(spec.city_in,:);
sharemat_10 = sharemat_10(:,[1 2 3 5]);
sharemat_07 = sharemat_07 ./ repmat(sum(sharemat_07,2),[1 5]);
sharemat_10 = sharemat_10 ./ repmat(sum(sharemat_10,2),[1 4]);

hhivec_07 = sum(sharemat_07.^2,2);
hhivec_10 = sum(sharemat_10.^2,2);

% Supermarkup vs. ABI share: 2007
scatter(sharemat_07(:,1),sm_07)
set(gca,'YLim',[0.0,2.5])
set(gca,'XLim',[0.1,0.8])
lbf = lsline;
lbf.Color = 'blue';
grid on
ylabel('')
xlabel('ABI Share','FontSize',20)
ylabel('Supermarkup','FontSize',20)
title('2007','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('smresw_abi_07'), 'pdf')
cd(path.code1)


% Supermarkup vs. Coors Share: 2007
scatter(sharemat_07(:,4),sm_07)
set(gca,'YLim',[0.0,2.5])
set(gca,'XLim',[0.0,0.3])
lbf = lsline;
lbf.Color = 'blue';
grid on
ylabel('')
xlabel('Coors Share','FontSize',20)
ylabel('Supermarkup','FontSize',20)
title('2007','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('smresw_coors_07'), 'pdf')
cd(path.code1)


% Supermarkup vs. ABI share: 2010
scatter(sharemat_10(:,1),sm_10)
set(gca,'YLim',[0.0,2.5])
set(gca,'XLim',[0.1,0.8])
lbf = lsline;
lbf.Color = 'blue';
grid on
ylabel('')
xlabel('ABI Share','FontSize',20)
ylabel('Supermarkup','FontSize',20)
title('2010','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('smresw_abi_10'), 'pdf')
cd(path.code1)


% Supermarkup vs. MillerCoors Share: 2010
scatter(sharemat_10(:,4),sm_10)
set(gca,'YLim',[0.0,2.5])
set(gca,'XLim',[0.2,0.8])
lbf = lsline;
lbf.Color = 'blue';
grid on
ylabel('')
xlabel('MillerCoors Share','FontSize',20)
ylabel('Supermarkup','FontSize',20)
title('2010','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('smresw_millercoors_10'), 'pdf')
cd(path.code1)

% Supermarkup vs. HHI: 2007
scatter(hhivec_07,sm_07)
set(gca,'YLim',[0.0,2.5])
set(gca,'XLim',[0.2,0.7])
lbf = lsline;
lbf.Color = 'blue';
grid on
ylabel('')
xlabel('HHI','FontSize',20)
ylabel('Supermarkup','FontSize',20)
title('2007','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('smresw_hhi_07'), 'pdf')
cd(path.code1)

% Supermarkup vs. HHI: 2010
scatter(hhivec_10,sm_10)
set(gca,'YLim',[0.0,2.5])
set(gca,'XLim',[0.2,0.7])
lbf = lsline;
lbf.Color = 'blue';
grid on
ylabel('')
xlabel('HHI','FontSize',20)
ylabel('Supermarkup','FontSize',20)
title('2010','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('smresw_hhi_10'), 'pdf')
cd(path.code1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure G4: Implied friction
%    - also assorted numbers cited in text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Load results at selected discount factor
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc pnash snash sm_10;
mc0=mc;
load sres_bind_2010 mc pnash snash smfinal gicc outdata
cd(path.code1)

firmidx = ids.firmid(ids.fiscid==2010);

% Specifications
lspec.skipcity = zeros(37,1);
lspec.fullvec = 2;
lspec.getder = 0;
lspec.scale = 100;
lspec.purpose = 'S';
lspec.icc_restrict = 0;
lspec.mcleader = 0;

idat.mc=mc0;
idat.pnash=mc;
idat.snash=mc;
idat.mc(ids.fiscid==2010)= mc;
idat.pnash(ids.fiscid==2010)=pnash;
idat.snash(ids.fiscid==2010)=snash;

smdevs = [];
mktvec = unique(ids.cdid(ids.fiscid==2010))';

[~,~,usedata,~,~,giccx] = f_loss_bind(sm_10,mktvec,vars,ids,idat,daugfile,spec,lspec,smdevs,[]);

pic = sum(usedata.picoll(firmidx==5));
pid = sum(usedata.pidev(firmidx==5));
pib = sum(usedata.pinash(firmidx==5));

etax = [0.26 (0.3:0.05:0.9)]';
fricx = zeros(size(etax));
for i=1:length(etax)
    fricx(i) = pic/(1-etax(i)) - pid - pib*etax(i)/(1-etax(i));
end
fricxrat = fricx ./ pic;

disp([etax fricxrat]);

%Writing to .txt file
cd(strcat(path.data2))
    TT = table([etax fricxrat]);
    writetable(TT,'impliedfrictions.txt','Delimiter',';');
cd(path.code1)

% Joint Identification
figure
plot(etax,fricxrat, 'k' )
    grid on
    ylabel('Implied Friction','FontSize',20)
    xlabel('Timing Parameter','FontSize',20)
    set(gcf,'color','white');
    set(gcf, 'PaperPosition', [0 0 5 5]); 
    set(gcf, 'PaperSize', [5 5]); 
cd(path.figr)
saveas(gcf, strcat('impliedfriction'), 'pdf')
cd(path.code1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


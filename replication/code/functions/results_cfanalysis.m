function results_cfanalysis(path)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function prepares Tables 6, 7, and 8, and Figures 6 and G3.
%
% Called by:
%   - main_v7.m
% Calls the following user-specified functions:
%   - main_data.m
%   - f_loss_bind.m
%   - f_inclusive.m
% Loads results created in the functions:
%   - combine_imputed.m
%   - results_cfmergers.m
%   - results_cfmmc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Baseline specification
dfx = 26;
df = dfx/100;
spec = main_spec(df);
spec.bysize = 0;  

% Grabbing data and defining fiscnum
[vars,ids] = main_data(path,spec);

firmid_11 = ids.firmid(ids.fiscid==2011);

firmid_10 = ids.firmid(ids.fiscid==2010);
firmid_10(ids.brndid(ids.fiscid==2010)==3) = 4;
firmid_10(ids.brndid(ids.fiscid==2010)==4) = 4;


% Demand results
cd(strcat(path.data2))
load daugfile daugfile;
cd(path.code1)
pcoefi = daugfile.pcoefi(ids.obsindemand==1,:);
deltanp = daugfile.deltanp(ids.obsindemand==1);
mu = daugfile.mu(ids.obsindemand==1,:);
ai = daugfile.ai(ids.obsindemand==1,:);
delta = daugfile.delta(ids.obsindemand==1);

% Results at selected discount factor
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc pnash snash sm_06 sm_07 sm_10 sm_11;
cd(path.code1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tables 6 and 8:
% Reading in merger counterfactual information and constructing table
%   - Changes in prices, m, profit, CS due to merger ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shell file 
resmat = zeros(19,5);


for col=1:5

    %Data for ABI/Modelo scenarios
    if col==1
        
        cd(strcat(path.data2))
        load cfscen_06 outdata sm_cf modsmfinal
        cd(path.code1)
        
    elseif col==2
        
        cd(strcat(path.data2))
        load cfscen_07 outdata sm_cf modsmfinal
        cd(path.code1)
        
    elseif col==3
        
        cd(strcat(path.data2))
        load cfscen_08 outdata sm_cf modsmfinal
        cd(path.code1)
                
    end
    
    % Filling in ABI/Modelo columns
    if col==1 || col==2 || col==3
    
        dprice = outdata.pcoll - vars.p_jt(ids.fiscid==2011);
        dBprice = outdata.pnash - pnash(ids.fiscid==2011);
        dshare = (outdata.scoll - vars.s_jt(ids.fiscid==2011))./ vars.s_jt(ids.fiscid==2011);

        %dcs_dpi = (sum(cs_new) - sum(cs_old)) / sum(dpi);
        
        % Supermarkups
        resmat(1,col) = mean(sm_cf);
        resmat(2,col) = mean(modsmfinal);
        
        % Change in supermarkups
        resmat(3,col) = mean(sm_cf) - mean(sm_11);
        resmat(4,col) = mean(modsmfinal);
        
        % Change in Bertrand price
        resmat(5,col) = mean(dBprice(firmid_11==1));
        resmat(6,col) = mean(dBprice(firmid_11==5));
        resmat(9,col) = mean(dBprice(firmid_11==2));
        
        % Change in total price
        resmat(10,col) = mean(dprice(firmid_11==1));
        resmat(11,col) = mean(dprice(firmid_11==5));
        resmat(14,col) = mean(dprice(firmid_11==2));
        
        % Change in share
        resmat(15,col) = 100*mean(dshare(firmid_11==1));
        resmat(16,col) = 100*mean(dshare(firmid_11==5));
        resmat(19,col) = 100*mean(dshare(firmid_11==2));
        
   
        
    end
    
    
    % Data for Miller/Coors columns is a bit more complicated--need all
    if col==4 || col==5 
         
        cd(strcat(path.data2))
        load cfscen_01 outdata sm_cf 
        orig_basedata = outdata;
        orig_sm = sm_cf;
        load cfscen_02 outdata sm_cf
        cd(path.code1)

        if col==4

            dprice = outdata.pcoll - orig_basedata.pcoll;
            dBprice = outdata.pnash - orig_basedata.pnash;
            dshare = (outdata.scoll - orig_basedata.scoll)./ orig_basedata.scoll;
            
            dsm = sm_cf - orig_sm;
            smlev = sm_cf;
            
        elseif col==5

            dprice = vars.p_jt(ids.fiscid==2010) - orig_basedata.pcoll;
            dBprice = pnash(ids.fiscid==2010) - orig_basedata.pnash;
            dshare = (vars.s_jt(ids.fiscid==2010) - orig_basedata.scoll)./ orig_basedata.scoll;
            
            dsm = sm_10 - orig_sm;
            smlev = sm_10;
            
        end
        



        % Post-merger supermarkups
        resmat(1,col) = mean(smlev);
        
        % Change in supermarkups
        resmat(3,col) = mean(dsm);
        
        % Change in Bertrand price
        resmat(5,col) = mean(dBprice(firmid_10==1));
        resmat(7,col) = mean(dBprice(firmid_10==5));
        resmat(8,col) = mean(dBprice(firmid_10==4));        
        resmat(9,col) = mean(dBprice(firmid_10==2));
        
        % Change in total price
        resmat(10,col) = mean(dprice(firmid_10==1));
        resmat(12,col) = mean(dprice(firmid_10==5));
        resmat(13,col) = mean(dprice(firmid_10==4));
        resmat(14,col) = mean(dprice(firmid_10==2));
        
        % Change in share
        resmat(15,col) = 100*mean(dshare(firmid_10==1));
        resmat(17,col) = 100*mean(dshare(firmid_10==5));
        resmat(18,col) = 100*mean(dshare(firmid_10==4));        
        resmat(19,col) = 100*mean(dshare(firmid_10==2));
        
       
    end
        
end

% Switching order of columns
resmatb = zeros(size(resmat));

resmatb(:,1) = resmat(:,5);
resmatb(:,2) = resmat(:,4);
resmatb(:,3:5) = resmat(:,1:3);


resmatbx = resmatb([3:14],:);

%Writing to .txt file
cd(strcat(path.figr))
    TT1 = table(round(resmatbx,2));
    writetable(TT1,'cfmerger.txt','Delimiter',';');
cd(path.code1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 7: Decomposition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specifications
lspec.skipcity = zeros(37,1);
lspec.fullvec = 2;
lspec.getder = 0;
lspec.scale = 100;
lspec.purpose = 'S';
lspec.icc_restrict = 0;
lspec.mcleader = 0;

% Timing parameter (ratio form)
etarat = 0.26/(1-0.26);

% Observed equilibrium in 2010 and 2011 (and earlier)
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc pnash snash sm_10 sm_11;
cd(path.code1)

% Simulated equilibrium with ABI/Modelo (no efficiencies) in 2011
cd(strcat(path.data2))
load cfscen_06 outdata modsmfinal sm_cf
outdata_sim_06 = outdata;
smfinal_sim_06 = sm_cf;
modsmfinal_sim_06 = modsmfinal;
cd(path.code1)

% Simulated equilibrium with ABI/Modelo (minor efficiencies) in 2011
cd(strcat(path.data2))
load cfscen_07 outdata modsmfinal sm_cf
outdata_sim_07 = outdata;
smfinal_sim_07 = sm_cf;
modsmfinal_sim_07 = modsmfinal;
cd(path.code1)

% Simulated equilibrium with ABI/Modelo (major efficiencies) in 2011
cd(strcat(path.data2))
load cfscen_08 outdata modsmfinal sm_cf
outdata_sim_08 = outdata;
smfinal_sim_08 = sm_cf;
modsmfinal_sim_08 = modsmfinal;
cd(path.code1)

% Simulated equilibrium without MillerCoors in 2010
cd(strcat(path.data2))
load cfscen_01 outdata sm_cf
outdata_sim_01 = outdata;
smfinal_sim_01 = sm_cf;
cd(path.code1)

% Simulated equilibrium with MillerCoors but no efficiencies in 2010
cd(strcat(path.data2))
load cfscen_02 outdata sm_cf
outdata_sim_02 = outdata;
smfinal_sim_02 = sm_cf;
cd(path.code1)

% Counterfactual of ABI/Modelo with coalition at pre-merger supermarkups
idx=ids;
idx.firmid(idx.brndid==5) = 1;
idx.firmid(idx.brndid==6) = 1;
idx.coalid = (idx.firmid==1|idx.firmid==5|idx.firmid==4);
  
smdevs = [];
idat.mc=mc;
idat.pnash=pnash;
idat.snash=snash;
idatx=idat;
idatx.pnash(ids.fiscid==2011) = outdata_sim_06.pnash;
idatx.snash(ids.fiscid==2011) = outdata_sim_06.snash;

idaty=idat;
idaty.mc(ids.fiscid==2011) = outdata_sim_07.mc;
idaty.pnash(ids.fiscid==2011) = outdata_sim_07.pnash;
idaty.snash(ids.fiscid==2011) = outdata_sim_07.snash;

idatz=idat;
idatz.mc(ids.fiscid==2011) = outdata_sim_08.mc;
idatz.pnash(ids.fiscid==2011) = outdata_sim_08.pnash;
idatz.snash(ids.fiscid==2011) = outdata_sim_08.snash;

mktvec = unique(ids.cdid(ids.fiscid==2011))';

specx=spec;
specx.modsm = 1;

%[~,~,decompdata1_abimod,~,~,~] = f_loss_bind(sm_11,mktvec,vars,idx,idatx,daugfile,specx,lspec,smdevs,zeros(size(sm_11)));

[~,~,decompdata0_abimod,~,~,~] = f_loss_bind(sm_11,mktvec,vars,ids,idat,daugfile,spec,lspec,[],[]);         %no merger
[~,~,decompdata1_abimod,~,~,~] = f_loss_bind(sm_11,mktvec,vars,idx,idatx,daugfile,specx,lspec,[],sm_11); %merger no efficiencies, old m (but apply to Modelo)
[~,~,decompdata2_abimod,~,~,~] = f_loss_bind(sm_11,mktvec,vars,idx,idaty,daugfile,specx,lspec,[],sm_11); %merger w/ minor efficiencies, old m (but apply to Modelo)
[~,~,decompdata3_abimod,~,~,~] = f_loss_bind(smfinal_sim_07,mktvec,vars,idx,idaty,daugfile,specx,lspec,[],modsmfinal_sim_07); %merger w/ minor efficiencies, new m
[~,~,decompdata4_abimod,~,~,~] = f_loss_bind(smfinal_sim_06,mktvec,vars,idx,idatx,daugfile,specx,lspec,[],modsmfinal_sim_06); %merger no efficiecies, new m
[~,~,decompdata5_abimod,~,~,~] = f_loss_bind(sm_11,mktvec,vars,idx,idatx,daugfile,spec,lspec,[],[]);      %merger no efficiecies, old m
[~,~,decompdata6_abimod,~,~,~] = f_loss_bind(smfinal_sim_08,mktvec,vars,idx,idatz,daugfile,specx,lspec,[],modsmfinal_sim_08); %merger major efficiecies, new m  

firmy_abimod = ids.firmid(ids.fiscid==2011);
firmx_abimod = idx.firmid(ids.fiscid==2011);


% Counterfactual of Miller/Coors with coalition at pre-merger supermarkups
idx=ids;
idx.firmid(idx.brndid==3) = 4;
idx.firmid(idx.brndid==4) = 4;
idx.coalid = (idx.firmid==1|idx.firmid==5|idx.firmid==4);

smdevs = [];

idat.mc=mc;
idat.pnash=mc;
idat.snash=mc;

idatx=idat;
idatx.mc(ids.fiscid==2010)=outdata_sim_01.mc;
idatx.pnash(ids.fiscid==2010) = outdata_sim_01.pnash;
idatx.snash(ids.fiscid==2010) = outdata_sim_01.snash;

idaty=idat;
idaty.mc(ids.fiscid==2010)=outdata_sim_02.mc;
idaty.pnash(ids.fiscid==2010) = outdata_sim_02.pnash;
idaty.snash(ids.fiscid==2010) = outdata_sim_02.snash;

idatz=idat;
idatz.mc(ids.fiscid==2010)= mc(ids.fiscid==2010);
idatz.pnash(ids.fiscid==2010)=pnash(ids.fiscid==2010);
idatz.snash(ids.fiscid==2010)=snash(ids.fiscid==2010);

mktvec = unique(ids.cdid(ids.fiscid==2010))';

specx=spec;


%[~,~,test,~,~,~] = f_loss_bind(sm_10,mktvec,vars,ids,idatz,daugfile,specx,lspec,smdevs,[]);


[~,~,decompdata0_mc,~,~,~] = f_loss_bind(smfinal_sim_01,mktvec,vars,idx,idatx,daugfile,spec,lspec,[],[]);  %no merger
[~,~,decompdata1_mc,~,~,~] = f_loss_bind(smfinal_sim_01,mktvec,vars,ids,idaty,daugfile,specx,lspec,[],[]); %merger no efficiencies; hold pre-merger m
[~,~,decompdata2_mc,~,~,~] = f_loss_bind(smfinal_sim_01,mktvec,vars,ids,idatz,daugfile,specx,lspec,[],[]); %merger w/ efficiencies; hold pre-merger m
[~,~,decompdata3_mc,~,~,~] = f_loss_bind(sm_10,mktvec,vars,ids,idatz,daugfile,specx,lspec,[],[]);          %merger w/ efficiencies; reequilibrated m
[~,~,decompdata4_mc,~,~,~] = f_loss_bind(smfinal_sim_02,mktvec,vars,ids,idaty,daugfile,specx,lspec,[],[]); %merger no efficiencies; reequilibrated m

firmy_mc = ids.firmid(ids.fiscid==2010);
firmx_mc = idx.firmid(ids.fiscid==2010);


% Output Table--decomposition
decomp2 = zeros(14,8);
for col = [1 2 3 4 5 6 7 8]
    
    if col==1
        usedata = decompdata0_mc;
        firmz = firmx_mc;
        bindingfirm = 4;
    elseif col==2
        usedata = decompdata1_mc;
        firmz = firmy_mc;
        bindingfirm = 5;
    elseif col==3
        usedata = decompdata2_mc;
        firmz = firmy_mc;
        bindingfirm = 5;
    elseif col==4
        usedata = decompdata3_mc;
        firmz = firmy_mc;
        bindingfirm = 5;
    elseif col==5
        usedata = decompdata0_abimod;
        firmz = firmy_abimod;
        bindingfirm = 5;
    elseif col==6
        usedata = decompdata1_abimod;
        firmz = firmx_abimod;
        bindingfirm = 5;
    elseif col==7
        usedata = decompdata2_abimod;
        firmz = firmx_abimod;
        bindingfirm = 5;  
    elseif col==8
        usedata = decompdata3_abimod;
        firmz = firmx_abimod;
        bindingfirm = 5;         
    end
    
    decomp2(1,col) = sum(usedata.picoll(firmz==1));
    decomp2(2,col) = sum(usedata.pidev(firmz==1));
    decomp2(3,col) = sum(usedata.pinash(firmz==1));
    decomp2(4,col) = etarat*sum(usedata.picoll(firmz==1));
    decomp2(5,col) = etarat*sum(usedata.pinash(firmz==1));
    decomp2(6,col) = (decomp2(4,col) - decomp2(5,col)) - (decomp2(2,col) - decomp2(1,col));
    
    decomp2(7,col) = sum(usedata.picoll(firmz==bindingfirm));
    decomp2(8,col) = sum(usedata.pidev(firmz==bindingfirm));
    decomp2(9,col) = sum(usedata.pinash(firmz==bindingfirm));
    decomp2(10,col) = etarat*sum(usedata.picoll(firmz==bindingfirm));
    decomp2(11,col) = etarat*sum(usedata.pinash(firmz==bindingfirm));
    decomp2(14,col) = (decomp2(10,col) - decomp2(11,col)) - (decomp2(8,col) - decomp2(7,col));

end

decomp2(12,:) = decomp2(7,:) +  decomp2(10,:);
decomp2(13,:) = decomp2(8,:) +  decomp2(11,:);

% Rescaling to millions
decomp2 = decomp2 / 1000000;

% Only tabled rows
decomp2 = decomp2(7:14,:);
decomp2 = decomp2([1 2 3 6 7 8],:);

%Writing to .txt file
cd(strcat(path.figr))
    TTD = table(round(decomp2,2));
    TTD2 = table(round(decomp2,3));
    writetable(TTD,'cfdecomp2.txt','Delimiter',';');
    writetable(TTD2,'cfdecomp2_threedigits.txt','Delimiter',';');
cd(path.code1)


% Output Table--WELFARE EFFECTS
weffects = zeros(7,8);
        
delta_base_10 = deltanp(ids.fiscid==2010)+daugfile.alpha*decompdata0_mc.pcoll;
delta_base_11 = deltanp(ids.fiscid==2011)+daugfile.alpha*decompdata0_abimod.pcoll;
mu_base_10 = mu(ids.fiscid==2010)-ai(ids.fiscid==2010,:).*repmat(vars.p_jt(ids.fiscid==2010),[1 size(ai,2)])+ai(ids.fiscid==2010,:).*repmat(decompdata0_mc.pcoll,[1 size(ai,2)]);
mu_base_11 = mu(ids.fiscid==2011,:);



cdid_10 = ids.cdid(ids.fiscid==2010);
temp1   = [ cdid_10(2:end); 0];
temp2   = (1:length(cdid_10))';
cdindex_10 = temp2(cdid_10~=temp1);

cdid_11 = ids.cdid(ids.fiscid==2011);
temp1   = [ cdid_11(2:end); 0];
temp2   = (1:length(cdid_11))';
cdindex_11 = temp2(cdid_11~=temp1);

cs_base_10 = f_inclusive(delta_base_10,mu_base_10,cdindex_10,pcoefi(ids.fiscid==2010),daugfile.rho,vars.msize(ids.fiscid==2010));
cs_base_11 = f_inclusive(delta_base_11,mu_base_11,cdindex_11,pcoefi(ids.fiscid==2011),daugfile.rho,vars.msize(ids.fiscid==2011));

%cs_old2 = f_inclusive(delta_base_11,mu_base_11,cdindex_11,pcoefi(ids.fiscid==2011),daugfile.rho,vars.msize(ids.fiscid==2011));


% Effect of Miller/Coors on Profit: Same m, no efficiencies vs. No Merger
weffects(1,1) = sum(decompdata1_mc.picoll) / sum(decompdata0_mc.picoll) -1;
weffects(2,1) = sum(decompdata1_mc.picoll(firmy_mc==1)) / sum(decompdata0_mc.picoll(firmy_mc==1)) -1;
weffects(4,1) = sum(decompdata1_mc.picoll(firmy_mc==5)) / sum(decompdata0_mc.picoll(firmy_mc==5)) -1;
weffects(5,1) = sum(decompdata1_mc.picoll(firmy_mc==2)) / sum(decompdata0_mc.picoll(firmy_mc==2)) -1;

delta_new_10 = deltanp(ids.fiscid==2010)+daugfile.alpha*decompdata1_mc.pcoll;
mu_new_10 = mu(ids.fiscid==2010)-ai(ids.fiscid==2010,:).*repmat(decompdata1_mc.pcoll,[1 size(ai,2)])+ai(ids.fiscid==2010,:).*repmat(decompdata1_mc.pcoll,[1 size(ai,2)]);
cs_new_10 = f_inclusive(delta_new_10,mu_new_10,cdindex_10,pcoefi(ids.fiscid==2010),daugfile.rho,vars.msize(ids.fiscid==2010));
weffects(7,1) = (sum(cs_new_10-cs_base_10)) / sum(decompdata1_mc.picoll - decompdata0_mc.picoll);

% Effect of Miller/Coors on Profit: Same m, efficiencies vs. No Merger
weffects(1,2) = sum(decompdata2_mc.picoll) / sum(decompdata0_mc.picoll) -1;
weffects(2,2) = sum(decompdata2_mc.picoll(firmy_mc==1)) / sum(decompdata0_mc.picoll(firmy_mc==1)) -1;
weffects(4,2) = sum(decompdata2_mc.picoll(firmy_mc==5)) / sum(decompdata0_mc.picoll(firmy_mc==5)) -1;
weffects(5,2) = sum(decompdata2_mc.picoll(firmy_mc==2)) / sum(decompdata0_mc.picoll(firmy_mc==2)) -1;

delta_new_10 = deltanp(ids.fiscid==2010)+daugfile.alpha*decompdata2_mc.pcoll;
mu_new_10 = mu(ids.fiscid==2010)-ai(ids.fiscid==2010,:).*repmat(decompdata2_mc.pcoll,[1 size(ai,2)])+ai(ids.fiscid==2010,:).*repmat(decompdata2_mc.pcoll,[1 size(ai,2)]);
cs_new_10 = f_inclusive(delta_new_10,mu_new_10,cdindex_10,pcoefi(ids.fiscid==2010),daugfile.rho,vars.msize(ids.fiscid==2010));
weffects(7,2) = (sum(cs_new_10-cs_base_10)) / sum(decompdata2_mc.picoll - decompdata0_mc.picoll);

% Effect of Miller/Coors on Profit: new m, no efficiencies vs. No Merger
weffects(1,3) = sum(decompdata4_mc.picoll) / sum(decompdata0_mc.picoll) -1;
weffects(2,3) = sum(decompdata4_mc.picoll(firmy_mc==1)) / sum(decompdata0_mc.picoll(firmy_mc==1)) -1;
weffects(4,3) = sum(decompdata4_mc.picoll(firmy_mc==5)) / sum(decompdata0_mc.picoll(firmy_mc==5)) -1;
weffects(5,3) = sum(decompdata4_mc.picoll(firmy_mc==2)) / sum(decompdata0_mc.picoll(firmy_mc==2)) -1;

delta_new_10 = deltanp(ids.fiscid==2010)+daugfile.alpha*decompdata4_mc.pcoll;
mu_new_10 = mu(ids.fiscid==2010)-ai(ids.fiscid==2010,:).*repmat(decompdata4_mc.pcoll,[1 size(ai,2)])+ai(ids.fiscid==2010,:).*repmat(decompdata4_mc.pcoll,[1 size(ai,2)]);
cs_new_10 = f_inclusive(delta_new_10,mu_new_10,cdindex_10,pcoefi(ids.fiscid==2010),daugfile.rho,vars.msize(ids.fiscid==2010));
weffects(7,3) = (sum(cs_new_10-cs_base_10)) / sum(decompdata4_mc.picoll - decompdata0_mc.picoll);

% Effect of Miller/Coors on Profit: Full vs. No Merger
weffects(1,4) = sum(decompdata3_mc.picoll) / sum(decompdata0_mc.picoll) -1;
weffects(2,4) = sum(decompdata3_mc.picoll(firmy_mc==1)) / sum(decompdata0_mc.picoll(firmy_mc==1)) -1;
weffects(4,4) = sum(decompdata3_mc.picoll(firmy_mc==5)) / sum(decompdata0_mc.picoll(firmy_mc==5)) -1;
weffects(5,4) = sum(decompdata3_mc.picoll(firmy_mc==2)) / sum(decompdata0_mc.picoll(firmy_mc==2)) -1;

delta_new_10 = deltanp(ids.fiscid==2010)+daugfile.alpha*decompdata3_mc.pcoll;
mu_new_10 = mu(ids.fiscid==2010)-ai(ids.fiscid==2010,:).*repmat(decompdata3_mc.pcoll,[1 size(ai,2)])+ai(ids.fiscid==2010,:).*repmat(decompdata3_mc.pcoll,[1 size(ai,2)]);
cs_new_10 = f_inclusive(delta_new_10,mu_new_10,cdindex_10,pcoefi(ids.fiscid==2010),daugfile.rho,vars.msize(ids.fiscid==2010));
weffects(7,4) = (sum(cs_new_10-cs_base_10)) / sum(decompdata3_mc.picoll - decompdata0_mc.picoll);

% Effect of ABI/Modelo on Profit: old m, no efficiencies 
weffects(1,5) = sum(decompdata5_abimod.picoll) / sum(decompdata0_abimod.picoll) -1;
weffects(3,5) = sum(decompdata5_abimod.picoll(firmy_abimod==1 | firmy_abimod==2)) / sum(decompdata0_abimod.picoll(firmy_abimod==1 | firmy_abimod==2)) -1;
weffects(4,5) = sum(decompdata5_abimod.picoll(firmy_abimod==5)) / sum(decompdata0_abimod.picoll(firmy_abimod==5)) -1;

delta_new_11 = deltanp(ids.fiscid==2011)+daugfile.alpha*decompdata5_abimod.pcoll;
mu_new_11 = mu(ids.fiscid==2011)-ai(ids.fiscid==2011,:).*repmat(decompdata5_abimod.pcoll,[1 size(ai,2)])+ai(ids.fiscid==2011,:).*repmat(decompdata5_abimod.pcoll,[1 size(ai,2)]);
cs_new_11 = f_inclusive(delta_new_11,mu_new_11,cdindex_11,pcoefi(ids.fiscid==2011),daugfile.rho,vars.msize(ids.fiscid==2011));
weffects(7,5) = (sum(cs_new_11-cs_base_11)) / sum(decompdata5_abimod.picoll - decompdata0_abimod.picoll);


% Effect of ABI/Modelo on Profit: new m, no efficiencies 
weffects(1,6) = sum(decompdata4_abimod.picoll) / sum(decompdata0_abimod.picoll) -1;
weffects(3,6) = sum(decompdata4_abimod.picoll(firmy_abimod==1 | firmy_abimod==2)) / sum(decompdata0_abimod.picoll(firmy_abimod==1 | firmy_abimod==2)) -1;
weffects(4,6) = sum(decompdata4_abimod.picoll(firmy_abimod==5)) / sum(decompdata0_abimod.picoll(firmy_abimod==5)) -1;

delta_new_11 = deltanp(ids.fiscid==2011)+daugfile.alpha*decompdata4_abimod.pcoll;
mu_new_11 = mu(ids.fiscid==2011)-ai(ids.fiscid==2011,:).*repmat(decompdata4_abimod.pcoll,[1 size(ai,2)])+ai(ids.fiscid==2011,:).*repmat(decompdata4_abimod.pcoll,[1 size(ai,2)]);
cs_new_11 = f_inclusive(delta_new_11,mu_new_11,cdindex_11,pcoefi(ids.fiscid==2011),daugfile.rho,vars.msize(ids.fiscid==2011));
weffects(7,6) = (sum(cs_new_11-cs_base_11)) / sum(decompdata4_abimod.picoll - decompdata0_abimod.picoll);

% Effect of ABI/Modelo on Profit: new m, minor efficiencies
weffects(1,7) = sum(decompdata3_abimod.picoll) / sum(decompdata0_abimod.picoll) -1;
weffects(3,7) = sum(decompdata3_abimod.picoll(firmy_abimod==1 | firmy_abimod==2)) / sum(decompdata0_abimod.picoll(firmy_abimod==1 | firmy_abimod==2)) -1;
weffects(4,7) = sum(decompdata3_abimod.picoll(firmy_abimod==5)) / sum(decompdata0_abimod.picoll(firmy_abimod==5)) -1;

delta_new_11 = deltanp(ids.fiscid==2011)+daugfile.alpha*decompdata3_abimod.pcoll;
mu_new_11 = mu(ids.fiscid==2011)-ai(ids.fiscid==2011,:).*repmat(decompdata3_abimod.pcoll,[1 size(ai,2)])+ai(ids.fiscid==2011,:).*repmat(decompdata3_abimod.pcoll,[1 size(ai,2)]);
cs_new_11 = f_inclusive(delta_new_11,mu_new_11,cdindex_11,pcoefi(ids.fiscid==2011),daugfile.rho,vars.msize(ids.fiscid==2011));
weffects(7,7) = (sum(cs_new_11-cs_base_11)) / sum(decompdata3_abimod.picoll - decompdata0_abimod.picoll);

% Effect of ABI/Modelo on Profit: new m, major efficiencies
weffects(1,8) = sum(decompdata6_abimod.picoll) / sum(decompdata0_abimod.picoll) -1;
weffects(3,8) = sum(decompdata6_abimod.picoll(firmy_abimod==1 | firmy_abimod==2)) / sum(decompdata0_abimod.picoll(firmy_abimod==1 | firmy_abimod==2)) -1;
weffects(4,8) = sum(decompdata6_abimod.picoll(firmy_abimod==5)) / sum(decompdata0_abimod.picoll(firmy_abimod==5)) -1;

delta_new_11 = deltanp(ids.fiscid==2011)+daugfile.alpha*decompdata6_abimod.pcoll;
mu_new_11 = mu(ids.fiscid==2011)-ai(ids.fiscid==2011,:).*repmat(decompdata6_abimod.pcoll,[1 size(ai,2)])+ai(ids.fiscid==2011,:).*repmat(decompdata6_abimod.pcoll,[1 size(ai,2)]);
cs_new_11 = f_inclusive(delta_new_11,mu_new_11,cdindex_11,pcoefi(ids.fiscid==2011),daugfile.rho,vars.msize(ids.fiscid==2011));
weffects(7,8) = (sum(cs_new_11-cs_base_11)) / sum(decompdata6_abimod.picoll - decompdata0_abimod.picoll);

weffects(1:6,:) = 100*weffects(1:6,:);

%Writing to .txt file
cd(strcat(path.figr))
    TTW = table(round(weffects,2));
    writetable(TTW,'cfmergerwelfare.txt','Delimiter',';');
cd(path.code1)

           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6: Multimarket competition  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Counterfactual with no MMC
cd(strcat(path.data2))
load cfscen_nommc sm ;
smx = sm;
sm_06x = smx(:,1);
sm_07x = smx(:,2);
sm_10x = smx(:,3);
sm_11x = smx(:,4);
cd(path.code1)

% Observed equilibrium
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc pnash snash sm_06 sm_07 sm_10 sm_11;
sm = [sm_06 sm_07 sm_10 sm_11];
cd(path.code1)

% Fleshing out information
lspec.skipcity = zeros(37,1);
lspec.fullvec = 2;
lspec.getder = 0;
lspec.scale = 100;
lspec.purpose = 'S';
lspec.icc_restrict = 0;
lspec.mcleader = 0;
idat.mc = mc;
idat.pnash = pnash;
idat.snash = snash;
mktvec_07 = unique(ids.cdid(ids.fiscid==2007))';
mktvec_10 = unique(ids.cdid(ids.fiscid==2010))';

firmid_07 = ids.firmid(ids.fiscid==2007);
firmid_10 = ids.firmid(ids.fiscid==2010);

[~,~,outdata_07x,~,~,~] = f_loss_bind(sm_07x,mktvec_07,vars,ids,idat,daugfile,spec,lspec,[],[]);
[~,~,outdata_07, ~,~,~] = f_loss_bind(sm_07, mktvec_07,vars,ids,idat,daugfile,spec,lspec,[],[]);
[~,~,outdata_10x,~,~,~] = f_loss_bind(sm_10x,mktvec_10,vars,ids,idat,daugfile,spec,lspec,[],[]);
[~,~,outdata_10, ~,~,~] = f_loss_bind(sm_10, mktvec_10,vars,ids,idat,daugfile,spec,lspec,[],[]);

% Average supermarkups 
[mean(sm) ; mean(smx) ; mean(sm-smx) ];

% Percentage effect on profit: 2007
picollx_a = sum(outdata_07x.picoll(firmid_07==1));
picollx_m = sum(outdata_07x.picoll(firmid_07==4));
picollx_c = sum(outdata_07x.picoll(firmid_07==5));
picoll_a  = sum(outdata_07.picoll(firmid_07==1));
picoll_m  = sum(outdata_07.picoll(firmid_07==4));
picoll_c  = sum(outdata_07.picoll(firmid_07==5));

[picoll_a picoll_m picoll_c ; picollx_a picollx_m picollx_c];

out1a = [picoll_a/picollx_a picoll_m/picollx_m picoll_c/picollx_c] -1;


% Percentage effect on profit: 2010
picollx_a = sum(outdata_10x.picoll(firmid_10==1));
picollx_m = sum(outdata_10x.picoll(firmid_10==5));
picoll_a  = sum(outdata_10.picoll(firmid_10==1));
picoll_m  = sum(outdata_10.picoll(firmid_10==5));

[picoll_a picoll_m  ; picollx_a picollx_m ];

out1b = [picoll_a/picollx_a picoll_m/picollx_m ] -1;

%Writing to .txt file
cd(strcat(path.figr))
    TTa = table(round(out1a,4));
    TTb = table(round(out1b,4));
    writetable(TTa,'mmc_profit_2007.txt','Delimiter',';');
    writetable(TTb,'mmc_profit_2010.txt','Delimiter',';');
cd(path.code1)


% Scatterplots of similarity: 2007
scatter(sm_07,sm_07x)
set(gca,'XLim',[0.0,2.5])
set(gca,'YLim',[0.0,2.5])
dline = refline(1,0);
dline.Color = 'black';
grid on
ylabel('')
xlabel('Pooled Slack Function','FontSize',20)
ylabel('Region-Specific Slack Function','FontSize',20)
title('Supermarkups in 2007','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('nommc_07'), 'pdf')
cd(path.code1)


% Scatterplots of similarity: 2010
scatter(sm_10,sm_10x)
set(gca,'XLim',[0.0,2.5])
set(gca,'YLim',[0.0,2.5])
dline = refline(1,0);
dline.Color = 'black';
grid on
ylabel('')
xlabel('Pooled Slack Function','FontSize',20)
ylabel('Region-Specific Slack Function','FontSize',20)
title('Supermarkups in 2010','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('nommc_10'), 'pdf')
cd(path.code1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure G3: MillerCoors as the leader
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Counterfactual with no MMC
cd(strcat(path.data2))
load cfscen_12 sm_cf ;
cd(path.code1)

% Observed equilibrium
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine sm_10 ;
sm = sm_10;
cd(path.code1)



% Scatterplots of similarity: 2010
scatter(sm_10,sm_cf)
set(gca,'XLim',[0.5,3.0])
set(gca,'YLim',[0.5,3.0])
dline = refline(1,0);
dline.Color = 'black';
grid on
ylabel('')
xlabel('Leader is ABI','FontSize',20)
ylabel('Leader is MillerCoors','FontSize',20)
%title('Supermarkups in 2010','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('mcleader_10'), 'pdf')
cd(path.code1)

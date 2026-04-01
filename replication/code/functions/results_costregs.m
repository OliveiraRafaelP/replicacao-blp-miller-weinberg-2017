function results_costregs(path)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function prepares runs cost regressions for the purposes of model 
% selection. It produces the following exhibits shown in the paper:
%   - Tables for discount factor selection: pooled ICs, non-pooled ICs, and
%     pooled-ICs with two supermarkups. (Tables 3, G1, G2)
%   - Scatterplot of 1-m vs. 2-m models (Figure G1)
%   - Scatterplot for pooled ICs vs. non-pooled ICs (Figure G2) 
%
% Called by:
%   - main_v7.m
% Calls the following user-specified functions:
%   - main_data.m
%   - main_spec.m
%   - f_loss_bind.m
%   - licols.m 
% Reads in data created in the following functions:
%   - impute_bertrand.m
%   - combine_imputed.m
%   - main_supply_bind_nopool.m
%   - main_supply_nonbind.m
% Regression coefficients are used in the following function:
%   - results_cfmergers.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLS regressions of implied costs on shifters
%   - Identifies plausible timing parameters 
%   - Output saved for creation of specification tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gettin the shifters
spec = main_spec(0);
[vars,ids] = main_data(path,spec);

% Cost vars & identifier for clustered SEs; w has city, region, time FEs
X = [vars.apost vars.w];
X = licols(X);
cv = ids.cityid;

%t = 1 : pooled 1-m model 
%t = 2 : pooled 2-m model
%t = 4 : non-pooled 1-m model

for t = [1 2 4]
    
   % Discount factors to loop over (0=Bertrand, 99=Non-Binding ICC)
    if t==1
        dfxin = [0 20 25 26 30 35 40 99];
    elseif t==2
        dfxin = [0 25 26 30 99];
    elseif t==3
        dfxin = [0 25 30 31 99];
    elseif t==4
        dfxin = [0 25 27 30 99];        
    end
    
    % Shell file and counter
    gammamat = zeros(4,length(dfxin));
    semat = gammamat;
    negmc = zeros(9,length(dfxin));
    counter = 1;
    
    % Looping through timing factors
    for dfx = dfxin
        
        % Setting specification
        cd(path.code1)
        df = dfx/100;
        spec = main_spec(df);
        
        if dfx==0
            
            cd(strcat(path.data2))
            load sres_bertrand mc;
            cd(path.code1)
            
        elseif dfx==99

            cd(strcat(path.data2))
            if (t==1 || t==4)
                load sres_nonbind mc sm_unc;
            elseif (t==2 || t==3)
                load sres_nonbind_x mc sm_unc;
            end
            cd(path.code1)
            
            sm_06 = reshape(sm_unc(:,1,:),size(sm_unc,1),size(sm_unc,3));
            sm_07 = reshape(sm_unc(:,2,:),size(sm_unc,1),size(sm_unc,3));
            sm_10 = reshape(sm_unc(:,3,:),size(sm_unc,1),size(sm_unc,3));
            sm_11 = reshape(sm_unc(:,4,:),size(sm_unc,1),size(sm_unc,3));
            
        else
            
            cd(strcat(path.data2,'/',spec.ssubfolder))
            if t==1
                load sres_bind_combine mc sm_06 sm_07 sm_10 sm_11;
            elseif t==2
                load sres_bind_x_combine mc sm_06 sm_07 sm_10 sm_11;
            elseif (t==4)
                load sres_bind_nopool mc sm;
                sm_06 = sm(:,1);
                sm_07 = sm(:,2);
                sm_10 = sm(:,3);
                sm_11 = sm(:,4);
            end
            
            cd(path.code1)
            
        end
        
       
        % OLS regression with clustered standard errors
        [gamma,se] = f_ols(mc,X,cv);
        
        
         %Saving cost parameters for baseline specification
        cd(strcat(path.data2))
        if t==1 && dfx==26
            save gammapar gamma se
        end
        cd(path.code1)
        
        
        
        %Filling in shell files
        gammamat(:,counter) = gamma(1:4);
        semat(:,counter) = se(1:4);
        
        if dfx~=0
            meansm = [mean(sm_06) mean(sm_07) mean(sm_10) mean(sm_11)];
            negmc(1:length(meansm),counter) = meansm';
        end
        negmc(end,counter) = length(mc(mc<0)) / length(mc);
        
        %Updating counter
        counter = counter+1;
    end
    
    TT1 = table(round(gammamat,3));
    TT2 = table(round(semat,3));
    temp1 = round(negmc,2);
    temp2 = round(negmc,3);
    temp1(end,:) = temp2(end,:);
    TT3 = table(temp1);    
    
    %Writing to .txt file
    cd(strcat(path.figr))
    if t==1
        writetable(TT1,'spec_df_coef_mod1.txt','Delimiter',';');
        writetable(TT2,'spec_df_serr_mod1.txt','Delimiter',';');
        writetable(TT3,'spec_df_mval_mod1.txt','Delimiter',';');
    elseif t==2
        writetable(TT1,'spec_df_coef_mod2.txt','Delimiter',';');
        writetable(TT2,'spec_df_serr_mod2.txt','Delimiter',';');
        writetable(TT3,'spec_df_mval_mod2.txt','Delimiter',';');
    elseif t==4
        writetable(TT1,'spec_df_coef_mod4.txt','Delimiter',';');
        writetable(TT2,'spec_df_serr_mod4.txt','Delimiter',';');
        writetable(TT3,'spec_df_mval_mod4.txt','Delimiter',';');
    end
    cd(path.code1)
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure G1 (scatterplots of 1-m and 2-m models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Getting the data
spec = main_spec(0);
[vars,ids] = main_data(path,spec);

% Pooled Specification w/ 2 m
cd(path.code1)
df = 26/100;
spec = main_spec(df);
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_x_combine mc sm_06 sm_07 sm_10 sm_11;
mc_p = mc;
sm_06_p = sm_06;
sm_07_p = sm_07;
sm_10_p = sm_10;
sm_11_p = sm_11;
cd(path.code1)

% Pooled Specification w/ 1 m
cd(path.code1)
df = 26/100;
spec = main_spec(df);
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc sm_06 sm_07 sm_10 sm_11;
mc_p1 = mc;
sm_06_p1 = sm_06;
sm_07_p1 = sm_07;
sm_10_p1 = sm_10;
sm_11_p1 = sm_11;
cd(path.code1)

% Supermarkup vectors
sm_p = [ sm_06_p ; sm_07_p ; sm_10_p ; sm_11_p];
sm_p1 = [ sm_06_p1 ; sm_07_p1 ; sm_10_p1 ; sm_11_p1];
sm_p_sma = sm_p(:,1);
sm_p_big = sm_p(:,2);

% Figure A
scatter(sm_p_sma,sm_p_big)
ylabel('24 Packs','FontSize',20)
xlabel('6 and 12 Packs','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]); 
set(gcf, 'PaperSize', [5 5]); 
xlim([0,3]); ylim([0,3]);
refline([1 0]);
grid on
cd(path.figr)
saveas(gcf, strcat('sm_scatter_bigsmall'), 'pdf')
cd(path.code1)

title('Panel A: Supermarkups','FontSize',20)
cd(path.figr)
saveas(gcf, strcat('sm_scatter_bigsmall_t'), 'pdf')
cd(path.code1)


% Figure B
scatter(sm_p1,mean([sm_p_sma sm_p_big],2))
ylabel('Two Supermarkups','FontSize',20)
xlabel('Single Supermarkup','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]); 
set(gcf, 'PaperSize', [5 5]); 
xlim([0,3]); ylim([0,3]);
refline([1 0]);
grid on
cd(path.figr)
saveas(gcf, strcat('sm_scatter_2v1'), 'pdf')
cd(path.code1)

title('Panel B: Supermarkups','FontSize',20)
cd(path.figr)
saveas(gcf, strcat('sm_scatter_2v1_t'), 'pdf')
cd(path.code1)



% Figure C
scatter(mc_p1,mc_p)
ylabel('Two Supermarkups','FontSize',20)
xlabel('Single Supermarkup','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]); 
set(gcf, 'PaperSize', [5 5]); 
%xlim([0,3]); ylim([0,3]);
refline([1 0]);
grid on
cd(path.figr)
saveas(gcf, strcat('mc_scatter_2v1'), 'pdf')
cd(path.code1)

title('Panel C: Marginal Costs','FontSize',20)
cd(path.figr)
saveas(gcf, strcat('mc_scatter_2v1_t'), 'pdf')
cd(path.code1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure G2 (scatterplots of pooled vs. non-pooled models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Getting the data
spec = main_spec(0);
[vars,ids] = main_data(path,spec);

% Pooled Specification w/ 1 m
cd(path.code1)
df = 26/100;
spec = main_spec(df);
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_combine mc sm_06 sm_07 sm_10 sm_11;
mc_p = mc;
sm_06_p = sm_06;
sm_07_p = sm_07;
sm_10_p = sm_10;
sm_11_p = sm_11;
cd(path.code1)


% Non-Pooled Specification
cd(path.code1)
df = 27/100;
spec = main_spec(df);
cd(strcat(path.data2,'/',spec.ssubfolder))
load sres_bind_nopool mc sm;
mc_np = mc(:);
mc_np = mc_np(mc_np~=0);
sm_06_np = sm(:,1);
sm_07_np = sm(:,2);
sm_10_np = sm(:,3);
sm_11_np = sm(:,4);
cd(path.code1)



%corr(mc_p,mc_np)


% Figure A
scatter(sm_07_p(:),sm_07_np(:))
ylabel('Regional IC Constraint','FontSize',20)
xlabel('Pooled IC Constraint','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]); 
set(gcf, 'PaperSize', [5 5]); 
xlim([0,3]); ylim([0,3]);
refline([1 0]);
grid on
cd(path.figr)
saveas(gcf, strcat('sm_07_scatter_pvnp'), 'pdf')
cd(path.code1)

title('Panel A: 2007 Supermarkups','FontSize',20)
cd(path.figr)
saveas(gcf, strcat('sm_07_scatter_pvnp_t'), 'pdf')
cd(path.code1)


% Figure B
scatter(sm_10_p(:),sm_10_np(:))
ylabel('Regional IC Constraint','FontSize',20)
xlabel('Pooled IC Constraint','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]); 
set(gcf, 'PaperSize', [5 5]); 
xlim([0,3]); ylim([0,3]);
refline([1 0]);
grid on
cd(path.figr)
saveas(gcf, strcat('sm_10_scatter_pvnp'), 'pdf')
cd(path.code1)

title('Panel B: 2010 Supermarkups','FontSize',20)
cd(path.figr)
saveas(gcf, strcat('sm_10_scatter_pvnp_t'), 'pdf')
cd(path.code1)


% Figure C
scatter(mc_p,mc_np)
ylabel('Regional IC Constraint','FontSize',20)
xlabel('Pooled IC Constraint','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]); 
set(gcf, 'PaperSize', [5 5]); 
%xlim([0,3]); ylim([0,3]);
refline([1 0]);
grid on
cd(path.figr)
saveas(gcf, strcat('mc_scatter_pvnp'), 'pdf')
cd(path.code1)

title('Panel C: Marginal Costs','FontSize',20)
cd(path.figr)
saveas(gcf, strcat('mc_scatter_pvnp_t'), 'pdf')
cd(path.code1)








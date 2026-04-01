function idgraph(path)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates Figure 2, an illustration of identification
%
% Called by:
%   - main_v7.m
% Calls the following user-specified functions:
%   - main_spec.m
%   - main_data.m
%   - f_loss_bind.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dfx = 26;
df = dfx/100;
spec = main_spec(df);
spec.bysize = 0;

%Loading demand results
cd(strcat(path.data2))
load daugfile daugfile
cd(strcat(path.code1))
    
% Main data 
[vars,ids] = main_data(path,spec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting simulation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smdevs0 = zeros(length(spec.city_in),1);

% Initial imputation options
lspec.getder = 0;       % get the numerical derivatives
lspec.fullvec = 0;      % evaluate holding fixed deviations (0) or based on full sm vector (1) 
lspec.icc_restrict = 1; % restrict ICC analysis to Coors/MillerCoors
lspec.nd2 = 1;          % 0=one-sided numerical derivatives (1=two-sided)
lspec.purpose = 'I';    % Imputing not simulating
lspec.skipcity = zeros(size(smdevs0,1),1);
lspec.mcleader = 0;

% Number of steps
steps=100;

% Shell files
idgraphdata_06  = zeros(length(ids.cdid(ids.fiscid==2006)),2,steps+1); 
idgraphdata_07  = zeros(length(ids.cdid(ids.fiscid==2007)),2,steps+1); 
idgraphdata_10  = zeros(length(ids.cdid(ids.fiscid==2010)),2,steps+1); 
idgraphdata_11  = zeros(length(ids.cdid(ids.fiscid==2011)),2,steps+1); 


% FY 2006!
disp('Starting 2006');
mktvec = unique(ids.cdid(ids.fiscid==2006))';
for i=0:steps
    disp(i);
    [~,~,od,~,~,~] = f_loss_bind(i/50,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs0);
    idgraphdata_06(:,:,i+1) = [od.pnash od.mc];
end
cd(strcat(path.data2))
save idgraphdata_06 idgraphdata_06 
cd(strcat(path.code1))

% FY 2007!
disp('Starting 2007');
mktvec = unique(ids.cdid(ids.fiscid==2007))';
for i=0:steps
    disp(i);
    [~,~,od,~,~,~] = f_loss_bind(i/50,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs0);
    idgraphdata_07(:,:,i+1) = [od.pnash od.mc];
end
cd(strcat(path.data2))
save idgraphdata_07 idgraphdata_07 
cd(strcat(path.code1))

% FY 2010!
disp('Starting 2010');
mktvec = unique(ids.cdid(ids.fiscid==2010))';
for i=0:steps
    disp(i);
    [~,~,od,~,~,~] = f_loss_bind(i/50,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs0);
    idgraphdata_10(:,:,i+1) = [od.pnash od.mc];
end
cd(strcat(path.data2))
save idgraphdata_10 idgraphdata_10 
cd(strcat(path.code1))

% FY 2011!
disp('Starting 2011');
mktvec = unique(ids.cdid(ids.fiscid==2011))';
for i=0:steps
    disp(i);
    [~,~,od,~,~,~] = f_loss_bind(i/50,mktvec,vars,ids,1,daugfile,spec,lspec,smdevs0);
    idgraphdata_11(:,:,i+1) = [od.pnash od.mc];
end
cd(strcat(path.data2))
save idgraphdata_11 idgraphdata_11 
cd(strcat(path.code1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firmid06 = ids.firmid(ids.fiscid==2006);
brndid06 = ids.brndid(ids.fiscid==2006);

firmid07 = ids.firmid(ids.fiscid==2007);
brndid07 = ids.brndid(ids.fiscid==2007);

firmid10 = ids.firmid(ids.fiscid==2010);
brndid10 = ids.brndid(ids.fiscid==2010);

firmid11 = ids.firmid(ids.fiscid==2011);
brndid11 = ids.brndid(ids.fiscid==2011);

sm = (0:steps)/50;

cd(strcat(path.data2))
load idgraphdata_06 idgraphdata_06 
load idgraphdata_07 idgraphdata_07 
load idgraphdata_10 idgraphdata_10 
load idgraphdata_11 idgraphdata_11 
cd(strcat(path.code1))

ysel = 7;

if ysel==6
    idgraphdata = idgraphdata_06;
    firmidx = firmid06;
    brndidx = brndid06;
elseif ysel==7
    idgraphdata = idgraphdata_07;
    firmidx = firmid07;
    brndidx = brndid07;
elseif ysel==10
    idgraphdata = idgraphdata_10;
    firmidx = firmid10;
    brndidx = brndid10;
elseif ysel==11
    idgraphdata = idgraphdata_11;
    firmidx = firmid11;
    brndidx = brndid11;
end


pnashy = permute(idgraphdata(:,1,:),[1 3 2]);
mcy = permute(idgraphdata(:,2,:),[1 3 2]);

pnashy_bud12 = mean(pnashy(brndidx==1,:));
pnashy_cor12 = mean(pnashy(brndidx==5,:));

mcy_bud12 = mean(mcy(brndidx==1,:));
mcy_cor12 = mean(mcy(brndidx==5,:));


% Identification Figures
figure
plot(sm,pnashy_bud12, 'k', ...
    sm,mcy_bud12, 'k-.')
grid on
legend('Bertrand Price','Marginal Cost', 'Location', 'southwest','FontSize',18)
ylabel('Dollars','FontSize',20)
xlabel('Supermarkup','FontSize',20)
%xlabel('')
title('Bud Light 12 Pack','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('idgraph_a_2007'), 'pdf')
cd(path.code1)

% Identification Figures
figure
plot(sm,pnashy_cor12, 'k', ...
    sm,mcy_cor12, 'k-.')
set(gca,'YLim',[9,16])
grid on
legend('Bertrand Price','Marginal Cost', 'Location', 'southwest','FontSize',18)
%ylabel('Price/Cost','FontSize',20)
xlabel('Supermarkup','FontSize',20)
%xlabel('')
title('Corona Extra 12 Pack','FontSize',20)
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
cd(path.figr)
saveas(gcf, strcat('idgraph_b_2007'), 'pdf')
cd(path.code1)

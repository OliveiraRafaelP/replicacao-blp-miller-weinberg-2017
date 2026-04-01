

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE SMALL SCANNER DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(path.data,'blp_beer_quarterly_brand_level'))    

% ID vectors 
id2 = data(:,1);
firmid = floor(id2/10000000000);
brndid = floor((id2-firmid*10000000000)/100000000);
sizeid = floor((id2-firmid*10000000000-brndid*100000000)/1000000);
cityid = floor((id2-firmid*10000000000-brndid*100000000-sizeid*1000000)/10000);
yearid = floor((id2-firmid*10000000000-brndid*100000000-sizeid*1000000-cityid*10000)/100);
montid = id2-firmid*10000000000-brndid*100000000-sizeid*1000000-cityid*10000-yearid*100;

% These observations fall outside the one-year period after the merger
obsin = (yearid<=3|yearid>=6)|(yearid==4&montid<=2)|(yearid==5&montid>=3);

% Defining the small sample for replication 
keeper = (sizeid==2 | sizeid==3); % No six-packs
keeper = keeper .* (brndid==1 | brndid==4 | brndid==5 | brndid==7 | brndid==13); % Main brands
keeper = keeper .* (cityid<6); % Five regions
keeper = keeper .*(1 - (brndid==5).*(sizeid==3)); % No Corona 24 pack
keeper = keeper .*(1 - (brndid==7).*(sizeid==3)); % No Heineken 24 pack.

% Need for demo draws below
yearcityid = grp2idx(yearid(keeper==1)*100+cityid(keeper==1));

% quantities increase bud by 1.5, miller by 1.5; perturb
quant = data(:,4) + (rand(length(keeper),1)-0.5);
quant(brndid==1) = 1.5*quant(brndid==1);
quant(brndid==13) = 1.5*quant(brndid==13);

% price; perturb
price = data(:,2) + 0.2*(rand(length(keeper),1)-0.5);
price(brndid==1) = 1.5*price(brndid==1);
price(brndid==13) = 1.5*price(brndid==13);

% output file
small_scanner=data(keeper==1,:);
small_scanner(:,2) = price(keeper==1);
small_scanner(:,4) = quant(keeper==1);

cd(strcat(path.data))
save small_scanner small_scanner;
cd(strcat(path.code))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE SMALL DEMO DRAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Selecting the draws 
demos = csvread(fullfile(path.data,'demosE.csv'),1,2);

% keep only the first five cities
cityindex = repmat([1:39]',[7 1]);
small_demosE = demos(cityindex<6,:);

%For use below, recalculating derivatives
demos2 = demos(yearcityid,:);

% Number of draws
ns= 500;

% Income, demeaned
inc  = demos2(:,1:ns);
inc  = inc - mean(mean(inc));

% Formatting observed demographics  
dfull = inc; 


save(fullfile(path.data,'small_demosE.mat'),'small_demosE');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE SMALL DEMO DRAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
% Loading demand parameters
cd(strcat(path.data,'/','RCNL2'))
load dres_gmm2
theta2=theta2_2;
alpha=theta1_2(1);
rho=rho_2;
theti = [1;2;3];
thetj = [2;2;2];
theta2w = full(sparse(theti,thetj,theta2));
clear k temp
cd(strcat(path.code))

% Recalculate derivatives with new, perturbed smalldata
cd(strcat(path.data))
load small_scanner
cd(strcat(path.code))

% ID vectors 
id2 = small_scanner(:,1);
firmid = floor(id2/10000000000);
brndid = floor((id2-firmid*10000000000)/100000000);
sizeid = floor((id2-firmid*10000000000-brndid*100000000)/1000000);
cityid = floor((id2-firmid*10000000000-brndid*100000000-sizeid*1000000)/10000);
yearid = floor((id2-firmid*10000000000-brndid*100000000-sizeid*1000000-cityid*10000)/100);
montid = id2-firmid*10000000000-brndid*100000000-sizeid*1000000-cityid*10000-yearid*100;

% price, quantity, shares
q_jt = small_scanner(:,4);
p_jt = small_scanner(:,2);
msize=small_scanner(:,7);
s_jt = q_jt./msize;

% x2 matrix
calor=small_scanner(:,6)/100;
calor = calor-mean(calor);
x2 = [p_jt ones(length(p_jt),1) calor];

% cdid
temp=num2str(id2);
[~,~,t1]=unique(str2num(temp(:,6:11)));
t2 = [0; t1(1:(length(t1)-1)) ];
cdid = ones(length(t1),1);
for i = 2:length(cdid)
    cdid(i) = cdid(i-1)+ 1*(t1(i)~=t2(i));
end
clear t1 t2 temp

% cdindex
temp1   = [ cdid(2:end); 0];
temp2   = (1:length(cdid))';
cdindex = temp2(cdid~=temp1);
clear temp1 temp2  

% getting mean consumer valuations, etc., with perturbed prices and shares
vars.p_jt = p_jt;
vars.s_jt = s_jt;
vars.ns = ns;
vars.dfull = dfull;
vars.x2 = x2;
temp = cumsum(s_jt);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
inshr = sum1(cdid,:);
vars.logcondshr = log(s_jt) - log(inshr);
vars.outshr = 1.0 - inshr;
daugfile.theta2w = theta2w;
daugfile.rho = rho;
idmatrix.cdindex = cdindex;
[delta,mu,ai] = rcnl_meanval(vars,idmatrix,daugfile);

% getting derivatives
derMat = zeros(8,8,140);
elasMat = zeros(8,8,140);
keepderivs = zeros(140,1);
for i = unique(cdid)'
   deltaM = delta(cdid==i,:);
   dfullM = dfull(cdid==i,:);
   x2M = x2(cdid==i,:);
   muM = mu(cdid==i,:);
   aiM = ai(cdid==i,:);
   pcoefiM = alpha+aiM;

   [sharei,scondi,sgroupi]=rcnl_indsh(exp(deltaM),exp(muM),rho,0,0,1);
   svec = mean(sharei,2);
   
   der1 = rcnl_der1(pcoefiM,sharei,scondi,sgroupi,rho);  
   temp1 = repmat(p_jt(cdid==i),[1 length(p_jt(cdid==i))])';
   temp2 = repmat(s_jt(cdid==i),[1 length(s_jt(cdid==i))]);
   elas = der1.*temp1./temp2;
   
   derMat(:,:,i) = der1;
   elasMat(:,:,i) = elas;
   
   keepderivs(i) = max(obsin(cdid==i));
   
   std(obsin(cdid==i));
   
   
end

% Need to drop the 3rd-dimension associated with first year after merger.
derMat_2 = derMat(:,:,keepderivs==1);
elasMat_2 = elasMat(:,:,keepderivs==1);


cd(strcat(path.data,'/','RCNL2'))
save small_dresgmm2 derMat_2 elasMat_2 theta2_2 theta1_2 rho_2
cd(strcat(path.code))


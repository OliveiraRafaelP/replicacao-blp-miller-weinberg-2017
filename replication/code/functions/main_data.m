function [vars,ids] = main_data(path,spec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function reads in the quarterly-level beer data used in Miller and
% Weinberg (2017) and creates the relevant data files for price leadership.
%
% Reads in the following datasets:
%   - demosE.csv
%   - blp_beer_quarterly_brand_level.mat
%
% Calls the following user-written functions:
%   - cr_dum
%
% Called by:
%   - f_daugment.m
%   - main_supply_bind.m
%   - main_supply_nonbind.m
%   - impute_bertrand.m
%   - main_supply_bind_nopool.m
%   - main_supply_bind_x.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINING RELEVANT OBSERVATIONS -- THIS DEFINES "OBSIN" VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load(fullfile(path.data,'blp_beer_quarterly_brand_level'))    
data2 = load(fullfile(path.data1,'small_scanner'));    
data = data2.small_scanner;

% ID vectors 
id2 = data(:,1);
firmid = floor(id2/10000000000);
brndid = floor((id2-firmid*10000000000)/100000000);
sizeid = floor((id2-firmid*10000000000-brndid*100000000)/1000000);
cityid = floor((id2-firmid*10000000000-brndid*100000000-sizeid*1000000)/10000);
yearid = floor((id2-firmid*10000000000-brndid*100000000-sizeid*1000000-cityid*10000)/100);
montid = id2-firmid*10000000000-brndid*100000000-sizeid*1000000-cityid*10000-yearid*100;

yearcityid = grp2idx(yearid*100+cityid);

% Exclude year immediately post-dating Miller/Coors.
if strcmp(spec.dfolder,'RCNL2')  || strcmp( spec.dfolder,'RCNL4')
    obsin = (yearid<=3|yearid>=6)|(yearid==4&montid<=2)|(yearid==5&montid>=3);
elseif strcmp(spec.dfolder,'RCNL1')  || strcmp( spec.dfolder,'RCNL3') || strcmp( spec.dfolder,'NL1')
    if min(yearid)==5
        yearid = yearid-4;
    end
    obsin = (yearid<=3|yearid>=6)|(yearid==4&montid<6)|(yearid==5&montid>=6);
end
obsintemp = obsin;


% Include only fiscal years 2006, 2007, 2010, 2011
if strcmp(spec.dfolder,'RCNL2')  || strcmp( spec.dfolder,'RCNL4') 
    fiscid = yearid + 1*(montid>=4);
elseif strcmp(spec.dfolder,'RCNL1')  || strcmp( spec.dfolder,'RCNL3')  || strcmp( spec.dfolder,'NL1')
    fiscid = yearid + 1*(montid>=10);
end
fiscid = fiscid + 2004;
obsin = obsin.*(ismember(fiscid,[2006 2007 2010 2011])==1);

% Include specified cities
if spec.includeallobs==0
    obsin = obsin.*(ismember(cityid,spec.city_in)==1);
elseif spec.includeallobs==1
    obsin = obsintemp;
end

% Use to select observation from saved demand data
obsindemand  = obsin(obsintemp==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADING MAIN DATASET WITH PRICES, SHARES, ID VECTORS, FIXED EFFECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load(fullfile(path.data,'blp_beer_quarterly_brand_level'))    
data2 = load(fullfile(path.data1,'small_scanner'));    
data = data2.small_scanner;

% ID vector 
id2 = data(obsin==1,1);
firmid = floor(id2/10000000000);
brndid = floor((id2-firmid*10000000000)/100000000);
sizeid = floor((id2-firmid*10000000000-brndid*100000000)/1000000);
cityid = floor((id2-firmid*10000000000-brndid*100000000-sizeid*1000000)/10000);
yearid = floor((id2-firmid*10000000000-brndid*100000000-sizeid*1000000-cityid*10000)/100);
montid = id2-firmid*10000000000-brndid*100000000-sizeid*1000000-cityid*10000-yearid*100;
prodid = grp2idx(brndid*100+sizeid);
dateid = grp2idx(yearid*100+montid);

%Consistency in yearid
if min(yearid)==5
    yearid = yearid-4;
end


% Fiscal year starts in October
if strcmp(spec.dfolder,'RCNL2')  || strcmp( spec.dfolder,'RCNL4') 
    fiscid = yearid + 1*(montid>=4);
elseif strcmp(spec.dfolder,'RCNL1')  || strcmp( spec.dfolder,'RCNL3')  || strcmp( spec.dfolder,'NL1')
    fiscid = yearid + 1*(montid>=10);
end
fiscid = fiscid + 2004;

%Fiscal year * city identifier
fisccity = grp2idx(cityid*100+fiscid);

%Coalition and leader identifiers
coalid = (firmid==1|firmid==5|firmid==4);
leadid = firmid==1;


% CDID is a consecutively numbered market identifier 
temp=num2str(id2);
[~,~,t1]=unique(str2num(temp(:,6:11)));
t2 = [0; t1(1:(length(t1)-1)) ];
cdid = ones(length(t1),1);
for i = 2:length(cdid)
    cdid(i) = cdid(i-1)+ 1*(t1(i)~=t2(i));
end
clear t1 t2 temp

% CDINDEX gives the last row of each market
temp1   = [ cdid(2:end); 0];
temp2   = (1:length(cdid))';
cdindex = temp2(cdid~=temp1);
clear temp1 temp2  

% Prices 
p_jt=data(obsin==1,2);

% Unit sales (144oz equivalents)
q_jt=data(obsin==1,4);

% Market size 
msize=data(obsin==1,7);

% Market shares
s_jt = q_jt./msize;

% Marketing categories
import  = ismember(firmid,2:3);

% Miles to brewery (not interacted with distance).
miles=data(obsin==1,5);
dist=data(obsin==1,8);
distbutfor=data(obsin==1,9);

% Calories
calor=data(obsin==1,6)/100;
calor = calor-mean(calor);

% Fixed effects matrices
brndfe = cr_dum(brndid);
prodfe = cr_dum(prodid);
yearfe = cr_dum(yearid);
montfe = cr_dum(montid);
cityfe = cr_dum(cityid);
datefe = cr_dum(dateid);
fiscfe = cr_dum(fiscid);
fcityfe= cr_dum(fisccity);


% Demand fixed effects specification -- yearfe comes last
fesd = [prodfe datefe(:,2:end)];

% Supply fixed effects specification 
fess = [prodfe cityfe(:,2:end) datefe(:,2:end)];

% Linear demand variables including fixed effects
x1 = sparse([p_jt fesd]) ;

% Demand variables receiving a random coefficient
if strcmp(spec.dfolder,'RCNL1')  || strcmp( spec.dfolder,'RCNL2')  || strcmp( spec.dfolder,'NL1')
    x2 = [p_jt ones(length(p_jt),1) calor];
elseif strcmp(spec.dfolder,'RCNL3')  || strcmp( spec.dfolder,'RCNL4') 
    x2 = [p_jt ones(length(p_jt),1) import calor sizeid];
end

%Marginal cost shifters
mcpost = (firmid==5) & yearid>=5;
mpost = ((brndid==11 | brndid==12 | brndid==13).*(yearid>=5)==1);
cpost = (brndid==3 | brndid==4).*(yearid>=5);
apost = (firmid==1) & yearid>=5;
w = [mpost cpost dist fess];

%Inside and outside good share
temp = cumsum(s_jt);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
inshr = sum1(cdid,:);
clear sum1 temp

% Placing in a structure
vars.p_jt = p_jt;
vars.s_jt = s_jt;
vars.logcondshr = log(s_jt) - log(inshr);
vars.msize = msize;
vars.inshr = inshr;
vars.outshr = 1.0 - inshr;
vars.logodds = log(s_jt) - log(vars.outshr);
vars.x1 = x1;
vars.x2 = x2;
vars.w = w;
vars.fesd = fesd;
vars.fess = fess;
vars.dist = dist;
vars.distbutfor = distbutfor;
vars.apost = apost;
vars.mpost = mpost;
vars.cpost = cpost;
vars.mcpost = mcpost;
vars.miles = miles;


% ID structure to simplify passing through to functions    
ids.cdid = cdid;
ids.firmid = firmid;
ids.cityid = cityid;
ids.yearid = yearid;
ids.brndid = brndid;
ids.prodid = prodid;
ids.montid = montid;
ids.sizeid = sizeid;
ids.cdindex = cdindex;
ids.dateid = dateid;
ids.fiscid = fiscid;
ids.fisccity = fisccity;
ids.coalid = coalid;
ids.leadid = leadid;
ids.obsindemand = obsindemand;

clear data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADING DATA WITH DEMOGRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns=500;

%Selecting the draws 
%demos = csvread(fullfile(path.data,'demosE.csv'),1,2);
demos2 = load(fullfile(path.data1,'small_demosE'));    
demos = demos2.small_demosE;

%Expanding the data, removing post-merger year. Sumdat3 used below.
demos2 = demos(yearcityid,:);
demos3 = demos2(obsin==1,:);

% Income, demeaned
inc  = demos3(:,1:ns);
inc  = inc - mean(mean(inc));

% Formatting observed demographics  
dfull = inc; 
demogr = dfull(cdindex,:);

% Placing in a structure
vars.dfull = dfull;
vars.demogr = demogr;
vars.ns = ns;

clear data 


  
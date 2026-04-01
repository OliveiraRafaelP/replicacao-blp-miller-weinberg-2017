function start_idcheck(path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This code conducts an exercise on identification using PLE with 
%%%% logit demand with (i) one region, and then (ii) multiple regions. The
%%%% steps include:
%%%%  (1) obtain structural parameters using a logit/Bertand calibration
%%%%  (2) simulate PLE as a constrained maximization problem
%%%%  (3) Take the PLE data and recover the supermarkup(s) and costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting seed
rng('default');

% Shell files
calmc=0;
recmc=0;
simm =0;
recm =0;
iccb =0;

% calibration based on share and margin (alt=0) or quality and cost (alt=1)
alt = 1;

% Starting timer
tic

%looping over scenarios
for i=1:100
    
    disp(i)
    
    % Number of (single product) firms
    nPlayer = 4 + floor(rand(1)*7);
        
    %conditional shares
    temp = rand([nPlayer 1]);
    sc = temp / sum(temp);
    sc = sort(sc,'descend');
    
    %outside good share with support over 0.2-0.5
    sout = 0.2 + rand()*0.3;
    
    %markup on first product with support over 0.25-0.75
    mark = 0.25 + rand()*0.5;
    
    %coalition structure
    coalid = [1 1 0 floor(2*rand(1,nPlayer-3)) ];
    
    %consider spectrum of discount factors
    delta = 0.2;
    while delta<0.85
        
        [mc1,mc2,smarkup1,smarkup2,bind] = idcheck(sc,sout,mark,coalid,delta,alt);
    
        calmc=[calmc;mc1];
        recmc=[recmc;mc2];
        simm =[simm;smarkup1];
        recm =[recm;smarkup2];
        iccb = [iccb;bind];
        
        delta = delta+0.1;

    end

end

timer = toc/60;

% Deleting first (shell) observation
calmc = calmc(2:end);
recmc = recmc(2:end);
simm = simm(2:end);
recm = recm(2:end);
iccb = iccb(2:end);

% Preparing data for saving
idcheckout.calmc = calmc;
idcheckout.recmc = recmc;
idcheckout.simm = simm;
idcheckout.recm = recm;
idcheckout.iccb = iccb;


% Analysis and .txt output
diffmc = calmc - recmc;
yx1 = abs(diffmc./recmc);
qtile1 = quantile(yx1,[0.50 0.75 0.90 0.95 0.975 0.99 1.00]);
test1 = 1 - sum( yx1 > 0.05 ) / length(diffmc);
test2 = 1 - sum( yx1 > 0.01 ) / length(diffmc);
test3 = 1 - sum( yx1 > 0.001 ) / length(diffmc);
bindfrac = sum(iccb) / length(iccb);

%Writing to .txt file
cd(strcat(path.figr))
TT = table([mean(recmc) qtile1 test1 test2 test3 bindfrac timer]);
writetable(TT,'idcheck_results.txt','Delimiter',';');
save idcheckout idcheckout
cd(strcat(path.code1))

function [calres]= cal_logit(sc,sout,m,p)

%unconditional shares
su = sc*(1-sout);

%price parameter
alpha = - 1 / m / (1-su(1));

% mean consumer valuation and mean quality
meanval = log(su) - log(sout);
xi = meanval - alpha*p;

%markups, marginal costs, profit
mark = - 1 ./ (1-su) / alpha;
mc = p - mark;
pi = mark.*su;

%output file
calres.alpha = alpha;
calres.xi = xi;
calres.mc = mc;
calres.mark = mark;
calres.pi = pi;
calres.p_n = p;
calres.s_n = su;

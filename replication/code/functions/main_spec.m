function [spec] = main_spec(df)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function provides some initial specification choices. Modifications
% are frequently made to specific elements.
%
% Called by:
%   - main_starter.m
%   - results_costregs.m
%   - results_baseanalysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spec.dfolder='RCNL2';   %Picks out NL1, RCNL1, RCNL2, RCNL3, RCNL4 from Miller/Weinberg 2017
spec.sfolder='v8smalldata';

spec.ssubfolder=sprintf('df_%d',df*1000);
%spec.ssubfolder=sprintf('df_%d',275);

% Discount factor
spec.df = df;

% Pool ICCs acroos regions?
spec.pool = 1;

% Supermarkup varies by size class? 
spec.bysize = 2;  % 1=full variation; 0=none; 2= (6 & 12) vs. 24/30

% Special Modelo supermarkup?
spec.modsm = 0;

% Markets to exclude and include --ADJUSTED FOR SMALL SAMPLE
% spec.city_out = [21 31]' ;
% spec.city_in = exclude( [1:39]',spec.city_out);
spec.city_out = []' ;
spec.city_in = exclude( [1:5]',spec.city_out);

% Include all observations?
spec.includeallobs = 0;


% Optimization options
spec.options.simp = optimset('GradObj','off','MaxIter',500,'MaxFunEvals',1000,...
    'Display','iter','TolFun',1e-4,'TolX',1e-4);

spec.options.cf = optimoptions(@fsolve,'Algorithm','levenberg-marquardt',...
    'MaxIter',500,'MaxFunEvals',2000,'Disp','off','StepTolerance',1e-12,'FunctionTolerance',1e-12);


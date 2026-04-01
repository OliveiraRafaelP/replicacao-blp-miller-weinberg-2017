function [pi]= leadpipinc(inc, p_n, coal, leadid, xi, mc, alpha, m, outside, prodid, optsol)

% Description: isolates the profit of only the leader firm when price rises
% by inc
% Inputs:
    % Let J be the number of products (excluding the outside good)
    % inc: scalar increase in price above Nash
    % p_n: (Jx1) vector of Nash equilibrium prices
    % coal: (number of firms x 1) dummy vector, equals 1 if that firm is in 
    % the coalition
    % leadid: (number of firms x1) dummy vector, equals 1 if that firm is
    % the leader
    % xi: (Jx1) vector of quality parameters
    % mc: (Jx1) vector of marginal costs
    % alpha: scalar price coefficient
    % m: scalar market size
    % outside: scalar normalization for outside good value
    % prodid: (Jx1) vector that identifies which firm goes with which
    % product
    % optsol: options for the fsolve routine
    
% Output:
    % pi: scalar profit for the lead firm (summed across all products) and
    % made negative

% Calculate profits for all coalition memebers
[~, pi_temp]= pipinc(inc, p_n, coal, xi, mc, alpha, m, outside, prodid, optsol);

% Isolate the leader's profit
pi = -pi_temp(leadid == 1);

function [zero]= nashfoc(price, pfix, pfixid, xi, mc, alpha, m, outside, prodid)

% Description: calculates the value of the nash logit FOC at a candidate
% price vector
% Inputs:
    % Let J be the number of products (excluding the outside good)
    % price: (J-Kx1) vector of prices, where K is the number of prices held
    % fixed
    % pfix: (Jx1) vector of prices where wish to hold some fixed (helpful 
    % for deviation calculations)
    % pfixid: (Jx1) dummy vector, equals 1 if prices are fixed
    % xi: (Jx1) vector of quality parameters
    % mc: (Jx1) vector of marginal costs
    % alpha: scalar price coefficient
    % m: scalar market size
    % outside: scalar normalization for outside good value
    % prodid: (Jx1) vector that identifies which firm goes with which
    % product
    
% Output:
    % zero: (Jx1) vector of FOC values
    
% Combine the fixed and non-fixed price vectors (pfixid = 1 means fixed)
p_temp = zeros(size(xi));
p_temp(pfixid == 1) = pfix(pfixid == 1);
p_temp(pfixid == 0) = price;

% Calculate share as an input into the FOC    
[share, ~, ~] = logit(p_temp, xi, mc, alpha, m, outside, prodid);

% Construct the margin/share derivative inputs to the FOC
deriv_temp = -alpha * (p_temp - mc) .* share;
deriv_temp = accumarray(prodid, deriv_temp);
deriv_temp = deriv_temp(prodid);

% Calculate the FOC
zero = deriv_temp(pfixid == 0) + alpha * (price - mc(pfixid == 0)) + 1;



function [p_c, pi]= pipinc(inc, p_n, coal, xi, mc, alpha, m, outside, prodid, optsol)

% Description: calculates the profits to a coalition for a given increase
% in prices (includes profits for firms not in the coalition)
% Inputs:
    % Let J be the number of products (excluding the outside good)
    % inc: scalar increase in price above Nash
    % p_n: (Jx1) vector of Nash equilibrium prices
    % coal: (number of firms x 1) dummy vector, equals 1 if that firm is in 
    % the coalition
    % xi: (Jx1) vector of quality parameters
    % mc: (Jx1) vector of marginal costs
    % alpha: scalar price coefficient
    % m: scalar market size
    % outside: scalar normalization for outside good value
    % prodid: (Jx1) vector that identifies which firm goes with which
    % product
    % optsol: options for the fsolve routine
    
% Output:
    % p_c: (Jx1) vector of prices, including both coalition and 
    % non-coalition firms
    % pi: (number of firms x 1) vector of profits, including both coalition 
    % and non-coalition firms

% Form the pfix vector 
pfix = p_n + inc;

% Deal with non-coalition firms
if min(coal) == 0
    % Calculate the prices for the non-coalition firms
    coal_temp = coal(prodid);
    f = @(p)nashfoc(p, pfix, coal_temp, xi, mc, alpha, m, outside, prodid);
    [p_fringe] = fsolve(f, mc(coal_temp == 0), optsol);

    % Combine coalition and non-coalition prices
    p_temp = zeros(size(xi));
    p_temp(coal_temp == 1) = pfix(coal_temp == 1);
    p_temp(coal_temp == 0) = p_fringe;
else
    p_temp = pfix;
end

% Calculate the profit for all firms
p_c = p_temp;
[~, ~, pi]= logit(p_c, xi, mc, alpha, m, outside, prodid);
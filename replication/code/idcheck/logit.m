function [share, pi, pif]= logit(price, xi, mc, alpha, m, outside, prodid)

% Description: calculates the logit demand shares and profits conditional
% on price
% Inputs:
    % Let J be the number of products (excluding the outside good)
    % price: (Jx1) vector of prices
    % xi: (Jx1) vector of quality parameters
    % mc: (Jx1) vector of marginal costs
    % alpha: scalar price coefficient
    % m: scalar market size
    % outside: scalar normalization for outside good value
    % prodid: (Jx1) vector that identifies which firm goes with which
    % product
    
% Output:
    % share: (Jx1) vector of product market shares
    % pi: (Jx1) vector of product-level profits
    % pif: (number of firms x 1) vector of firm-level profits

% Calculate Market Shares    
num = exp(xi + alpha * price);
denom = exp(outside) + sum(num);

share = num./denom;

% Calculate product-level profit
pi = ((price - mc) .* share) * m;

%Calculate firm-level profit
pif = accumarray(prodid, pi);
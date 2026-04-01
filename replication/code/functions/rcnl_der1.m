function der = rcnl_der1(pcoefi,sharei,scondi,sgroupi,rho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function obtains a matrix that contains the first derivatives
% of market share with respect to price.  It integrates numerically over
% the derivatives of the choice probabilities of the individuals
% considered.  Speed is optimized through the use of 3D matrix
% manipulation.
%
% Calculations are done for market cdid==m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of products in the market
J = size(pcoefi,1);

%Consumer-specific choice probabilities and price coefficients. 
su = sharei;
sc = scondi;
sg = sgroupi;
ai  = pcoefi;

%Rearranging dimensionality of choice probabilities 
su1 = permute(su,[1 3 2]);
adj = 1+rho/(1-rho)*(1./sg);

%Calculating derivatives
own   = (1/(1-rho))*mean(ai.*su.*(1-(1-rho)*su-rho*sc),2);
cross = -mean(repmat(permute(ai.*adj,[1 3 2]),[1 J 1]).*mtimesx(su1,su1,'T'),3);
der = cross - diag(diag(cross)) + diag(own);


 

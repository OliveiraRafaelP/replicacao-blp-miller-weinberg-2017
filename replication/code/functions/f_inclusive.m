
function [cs,csi] = f_inclusive(delta,mu,cdindex,pcoefi,rho,msize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtains the average inclusive value in each market, converted to be in
% dollars units, and scaled by market size.  Thus, comparing two scenarios, 
% the difference in the output is the compensating variation of Small and
% Rosen (1981), and reflects changes in consumer surplus.  Inputs include 
% the mean valuation (delta), the consumer-specific deviations (mu), the
% idmatrix from which we obtain cdid, cdindex, and the matrix of
% consumer-specific price coefficients pcoefi.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delta=delta_1;
% mu=mu_1;
% alphai=ai_1;
% pcoefi=pcoefi_1;
% cdindex=idmatrix.cdindex;
% msize=vars.msizelong;
% rho = rho_1;


v = exp((repmat(delta,[1 size(mu,2)])+mu)/(1-rho));
temp = cumsum(v);
Dg = temp(cdindex,:);
Dg(2:size(Dg,1),:) = diff(Dg);
D = exp(0) + Dg.^(1-rho);
alphai   = pcoefi(cdindex,:);
csi = (1./abs(alphai)).*log(D);
cs = mean(csi,2).*msize(cdindex); 


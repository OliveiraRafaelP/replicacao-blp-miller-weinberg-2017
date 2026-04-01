function [su,sc,sg] = rcnl_indsh(expmval,expmu,rho,cdid,cdindex,loopflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the choice probabilities of each individual and
% returns them in matrix form.  The model is RCNL, which collapses to
% standard RC logit if rho=0.  The choice probabilities returned are
% unconditional (su), conditional on the inside good (sc), and the inside
% group choice probability (sg).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns = size(expmu,2);

% For use looping over individual markets (cdid).
if loopflag==1
    
    eg = (expmu.*kron(ones(1,ns),expmval)).^(1/(1-rho));
    temp = cumsum(eg);
    sum1 = temp(end,:);
    denom1 = 1./(sum1);             %logit share denominator (no outside good)
    incValIn1 = (1-rho).*log(sum1);   %inclusive value of inside goods
    incValAll1 = log(1+exp(incValIn1)); %inclusive value of all goods
    denom = repmat(denom1,[size(eg,1) 1]);
    incValIn = repmat(incValIn1,[size(eg,1) 1]);
    incValAll = repmat(incValAll1,[size(eg,1) 1]);
    su = eg.*exp(incValIn)./exp(incValIn./(1-rho))./exp(incValAll); %unconditional choice probabilities
    sg = exp(incValIn)./exp(incValAll); %inside good choice probabilities
    sc = eg.*denom;                     %conditional choice probabilities
    
% No Looping
elseif loopflag==0

    eg = (expmu.*kron(ones(1,ns),expmval)).^(1/(1-rho));
    temp = cumsum(eg);
    sum1 = temp(cdindex,:);
    sum1(2:size(sum1,1),:) = diff(sum1);
    denom1 = 1./(sum1);               %logit share denominator (no outside good)
    incValIn1 = (1-rho).*log(sum1);   %inclusive value of inside goods
    incValAll1 = log(1+exp(incValIn1)); %inclusive value of all goods
    denom = denom1(cdid,:);
    incValIn = incValIn1(cdid,:);
    incValAll = incValAll1(cdid,:);    
    su = eg.*exp(incValIn)./exp(incValIn./(1-rho))./exp(incValAll); %unconditional choice probabilities
    sg = exp(incValIn)./exp(incValAll); %inside good choice probabilities
    sc = eg.*denom;                     %conditional choice probabilities
    
end




function f = f_ownMat(vec)

%###################################################################
% Function converts firm assignment into ownership matrix
%   Also can be used to dummy out observations from different markets
%   Inputs include:
%     vec= vector assigning products to firms (or markes)
%   Output includes:
%     (unnamed) ownership/market matrix

C=repmat(vec,1,size(vec,1));

%f = sparse(C==C');
f = C==C';



% % Test data
% vec=  [ 1 1 2 2 4 ]'
% 
% % This rejected version is slower
% y  = vec*vec';
% z1 = (ones(length(vec),1)*(vec')) .* (ones(length(vec),1)*(vec'));
% f  = (y==z1) ;





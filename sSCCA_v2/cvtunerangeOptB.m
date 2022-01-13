%algorithm for range of tuning parameters

function [Tauxrange, Tauyrange]= cvtunerangeOptB(X,Y,ncancorr,method)
%%%%Inputs
%X:         n x p dataset with network information, where rows are samples, and
%           columns are variables
%Y:         n x q dataset without network information, where rows are samples,
%           and columns are variables
%ncancorr:  number of desired canonical correlation vectors
%method:    Fused or Grouped
%%%%%Outputs
%TauXrange: a vector with tuning range parameters for X
%TauYrange: a vector with tuning range parameters for Y
%obtain first nonsparse input
[tildeA, tildeB, ~,Sig11invSig12,Sig22invSig21]=mynonsccaOptB(X,Y);

[n,p]=size(X);
[n,q]=size(Y);
mymaxlambdaalpha=NaN(ncancorr,1);
mymaxlambdabeta=NaN(ncancorr,1);
for j=1:ncancorr
   mymaxlambdaalpha(j,1)= norm(Sig11invSig12*tildeB(:,j),Inf);
   mymaxlambdabeta(j,1) = norm(Sig22invSig21*tildeA(:,j),Inf);
   if strcmp(method,'Fused')
   upperboundA=mymaxlambdaalpha(j,1);
   upperboundB= mymaxlambdabeta(j,1);
   Tauxrange(j,:)=[upperboundA*sqrt(log(p)/n)  upperboundA];
   Tauyrange(j,:)=[upperboundB*sqrt(log(q)/n) upperboundB];
   elseif strcmp(method,'Grouped')
   upperboundA=mymaxlambdaalpha(j,1)/1.6; 
   upperboundB= mymaxlambdabeta(j,1)/1.4; 
   Tauxrange(j,:)=[upperboundA*sqrt(log(p)/n)  upperboundA];
   Tauyrange(j,:)=[upperboundB*sqrt(log(q)/n) upperboundB]; 
   end
end

%Sandra E. Safo
%All rights reserved
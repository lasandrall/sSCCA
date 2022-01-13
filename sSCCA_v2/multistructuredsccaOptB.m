function [maxcorr,hatalpha,hatbeta]= multistructuredsccaOptB(X,Y,ncancorr,TauX,TauY,edgesX,edgesY,weightsX,weightsY,method,mygamma,myeta,asnormalize);

%--------------------------------------------------------------------------
%multistructuredsccaOptB.m: function to estimate structured sparse CCA for 
%given tuning parameters using Option B in published manuscript. 
%Use this option if p/q is large (say p/q>= 500)
%--------------------------------------------------------------------------

%If you want to select optimal tuning parameters, use
%cvmultistructuredsccaOptB.m

%USAGE
%[maxcorr,hatalpha,hatbeta]= multistructuredsccaOptB(X,Y,ncancorr,TauX,TauY,edgesX,edgesY,weightsX,weightsY,method,mygamma,myeta,asnormalize);
%[maxcorr,hatalpha,hatbeta]= multistructuredsccaOptB(X,Y,ncancorr,TauX,TauY,edgesX,edgesY,weightsX,weightsY,method);%defualt settings
%see test_structuredsccaOptB.m for examples

% This function depends on the CVX package. Please download at
%http://cvxr.com/cvx/download/


%%%%Inputs
%X:         n x p dataset with network information, where rows are samples, and
%           columns are variables
%Y:         n x q dataset without network information, where rows are samples,
%           and columns are variables
%ncancorr:  number of desired canonical correlation vectors
%TauX:      ncancorr*iteration by 3 matrix. [iteration ncancorr TauX];
%TauY:      ncancorr*iteration by 3 matrix. [iteration ncancorr TauY];           
%edgesX:    N x 2 matrix giving the pairwise connection of the variables in
%           X
%edgesX:    M x 2 matrix giving the pairwise connection of the variables in
%           Y
%weightsX:  p x 1 vector describing the degree of each variable in X. No zero
%           allowed. weights =d +1, where d is the degree for each variable.
%weightsY:  q x 1 vector describing the degree of each variable in Y. No zero
%           allowed. weights =d +1, where d is the degree for each variable.
%method:    Fused or Grouped
%mygamma:   gamma for Grouped. If method is Fused, this is not used. 
%           default is 2;
%myeta:       balances selection of singletons and networks. 0<myeta <1.
%           default is 0.5
%asnormalize: True or False. If True, data will be normalized to have mean 0
%                   and variance one for each variable. Default is True.

%output
%hatalpha:  estimated structured sparse CCA vector(s) for X
%hatbeta:   estimated structured sparse CCA vector(s) for Y
%maxcorr:   estimated canonical correlation coefficient(s)


%set defaults;
if(nargin <11)
    mygamma=2;
    myeta=0.5;
    asnormalize='True';
end

if(nargin <12)
    myeta=0.5;
    asnormalize='True';
end

if(nargin <13)
   asnormalize='True';
end

nX=size(X,1);
nY=size(Y,1);

if (nX~=nY)
    error('X and Y have different number of observations');
end


if strcmp(asnormalize,'True')
   Xold=mynormalize(X);
   Yold=mynormalize(Y);
elseif(strcmp(asnormalize,'False'))
   Xold=X;
   Yold=Y;   
end

X=Xold;
Y=Yold;

%prepare inputs
[tildealpha, tildebeta, tilderho,~,~ ]=mynonsccaOptB(Xold,Yold);
%initialize
tol=10^(-4);
mybeta1=tildebeta(:,1:ncancorr);
myalpha1=tildealpha(:,1:ncancorr);
myrho1=tilderho(1:ncancorr);
hatalpha=[]; hatbeta=[]; maxcorr=[];
for j=1:ncancorr
    if(max(TauX(TauX(:,2)==j,1))==1)
        maxiteration=5;
        TauXj=repmat(TauX(TauX(:,2)==j,[1,3]),maxiteration,1);   
        TauXj(:,1)= reshape(repmat(1:maxiteration,1,1),[1*maxiteration,1]);
        TauYj=repmat(TauY(TauY(:,2)==j,[1,3]),maxiteration,1);   
        TauYj(:,1)= reshape(repmat(1:maxiteration,1,1),[1*maxiteration,1]);
    else
        TauXj=TauX(TauX(:,2)==j,[1,3]);
        TauYj=TauY(TauY(:,2)==j,[1,3]);
        maxiteration=max(TauXj(:,1));
    end
    iter=0;
    diffalpha=1;
    diffbeta=1;
  if(j>1)
        ProjmX=Xold*(hatalpha(:,1:j-1)/(hatalpha(:,1:j-1)'* hatalpha(:,1:j-1)))* hatalpha(:,1:j-1)';  
        ProjmX(isnan(ProjmX))=0;
        X=Xold-ProjmX;  
        ProjmY=Yold*(hatbeta(:,1:j-1)/(hatbeta(:,1:j-1)'* hatbeta(:,1:j-1))* hatbeta(:,1:j-1)'); 
        ProjmY(isnan(ProjmY))=0;
        Y=Yold-ProjmY;    
  end

while ( ((iter<maxiteration) && and( (diffalpha > tol) , (diffbeta >tol)) ) )
       
   iter=iter+1;       
   if(iter==1)
   mybetaold=mybeta1(:,j);
   myalphaold=myalpha1(:,j);
   tilderhoold=myrho1(j);
   else
   mybetaold=mybeta(:,1);
   myalphaold=myalpha(:,1);
   tilderhoold=myrho(1);
   end
   [myalpha, mybeta,myrho]= structuredsccaOptB(X,Y,mybetaold, myalphaold,...
   tilderhoold,TauXj(TauXj(:,1)==iter,2),TauYj(TauYj(:,1)==iter,2),edgesX,edgesY,weightsX,weightsY,mygamma,myeta,method);          

        %convergence
        if( or( min(sum(abs(myalpha)))==0, min(sum(abs(mybeta)))==0))
             myalpha=myalphaold;
             mybeta=mybetaold;          
            break
        end    
        
        if (abs(tilderhoold-myrho) <= 10^-3)
                 myalpha=myalphaold;
                 mybeta=mybetaold;
                 break
        end
        
        
        diffalpha=norm(abs(myalpha) - abs(myalphaold),Inf);
        diffbeta=norm(abs(mybeta) -   abs(mybetaold),Inf);

end


maxcorr(j,1)=abs(corr(Xold*myalpha,Yold*mybeta));
hatalpha(:,j)=myalpha;
hatbeta(:,j)=mybeta;
end

%display some results
NonZeroAlpha=sum(hatalpha~=0);
NonZeroBeta=sum(hatbeta~=0);
dNZalpha=['Num non-zeros alpha: ', num2str(NonZeroAlpha)];
disp(dNZalpha);
dNZbeta=['Num non-zeros beta: ', num2str(NonZeroBeta)];
disp(dNZbeta);
dcorr=['Corr(X\alpha,Y\beta): ', num2str(maxcorr')];
disp(dcorr);
d=['Structured Sparse CCA Method used is: ', method, ' Option B'];
disp(d);


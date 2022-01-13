function cvout=cvmultistructuredsccaOptA(X,Y,ncancorr,edgesX,edgesY,weightsX,weightsY,method,mygamma,myeta,nfolds,ngrid,asnormalize);

%--------------------------------------------------------------------------
%cvmultistructuredsccaOptA.m: cross validation approach to
%select optimal tuning parameters for structured sparse CCA with Option 
%B in published manuscript. 
%Use this option if p/q is small to moderate (say p/q< 500)
%--------------------------------------------------------------------------
%
%DESCRIPTION:
%Function performs nfolds cross validation to select
%optimal tuning parameter, which is then applied on whole data to obtain
%canonical correlation vectors. 
%If you want to apply optimal tuning parameters to testing data, use
%multistructuredsccaOptA.m

%USAGE
%cvout=cvmultistructuredsccaOptA(X,Y,ncancorr,edgesX,edgesY,weightsX,weightsY,method,mygamma,myeta,nfolds,ngrid,normalize);
%cvout=cvmultistructuredsccaOptA(X,Y,ncancorr,edgesX,edgesY,weightsX,weightsY,method);%uses
%defualt settings
%see test_structuredsccaOptA.m for examples


%DATE: APRIL 14,2017
%
%AUTHORS
%MATLAB CODE WAS WRITTEN BY SANDRA E. SAFO (seaddosafo@gmail.com)
%
%REFERENCES
%Integrative analysis of transcriptomic and metabolomic
% data via sparse canonical correlation analysis with incorporation 
%of biological information (2017). Sandra E. Safo, Shuzhao Li, and Qi Long
% To appear in Biometrics


% This function depends on the CVX package. Please download at
%http://cvxr.com/cvx/download/


%%%%Inputs
%X:         n x p dataset with network information, where rows are samples, and
%           columns are variables
%Y:         n x q dataset without network information, where rows are samples,
%           and columns are variables
%ncancorr:  number of desired canonical correlation vectors
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
%nfolds:    number of cross validation folds. default is 5.
%ngrid:     number of tuning parameters for X and Y. Default is 10 ;
%asnormalize: True or False. If True, data will be normalized to have mean 0
%                   and variance one for each variable. Default is True.

%tuning parameter selected by cross searching
%first fix Tauy and search accross Taux; select the best tuning parameter
%next with the selected tuning parameter Taux, search across Tauy and
%choose the best. If multiple Taux or Tauy are optimal, last one is chosen.
%optimal tuning chosen using | |mean CVTrain| - |mean CVTest| |

%output
%cvout.optTauX:   optimal tuning parameter for X
%cvout.optTauY:   optimal tuning parameter for Y
%cvout.hatalpha:  estimated structured sparse CCA vector for X
%cvout.hatbeta:   estimated structured sparse CCA vector for Y
%cvout.maxcorr:   estimated canonical correlation coefficient
%cvout.method:    method used, Fused or Grouped
%cvout.tunerange: tuning range used


%set defaults;
if(nargin <9)
    mygamma=2;
    myeta=0.5;
    nfolds=5;
    ngrid=10;
    asnormalize='True';
end

if(nargin <10)
    myeta=0.5;
    nfolds=5;
    ngrid=10;
    asnormalize='True';
end

if(nargin <11)
    nfolds=5;
    ngrid=10;
    asnormalize='True';
end

if(nargin <12)
    ngrid=10;
    asnormalize='True';
end

if(nargin <13)
    asnormalize='True';
end

nX=size(X,1);
nY=size(Y,1);
n=nX;
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

tol=10^(-4);
maxiteration=5;

X=Xold;
Y=Yold;

%obtain tuning parameter range for X and Y;
[Tauxrange, Tauyrange]= cvtunerangeOptA(X,Y,ncancorr,method);


myoptTau=[];Tauxvec=[]; Tauyvec=[];

foldid = randsample([repmat(1:nfolds,1,floor(n/nfolds)) 1:mod(n,nfolds)],n);

for j=1:ncancorr;
        iter=0;
        diff_alpha=1;
        diff_beta=1;
        
        Tauxvec=2.^(linspace(log(Tauxrange(j,1))./log(2), log(Tauxrange(j,2))./log(2), ngrid+1));
        Tauyvec=2.^(linspace(log(Tauyrange(j,1))./log(2), log(Tauyrange(j,2))./log(2), ngrid+1));
        Tauxvec=Tauxvec(1:ngrid);
        Tauyvec=Tauyvec(1:ngrid);
        if(j>1)
            ProjmX=Xold*(myalphamat(:,1:j-1)/(myalphamat(:,1:j-1)'* myalphamat(:,1:j-1)))* myalphamat(:,1:j-1)';  
            ProjmX(isnan(ProjmX))=0;
            X=Xold-ProjmX;  
    %         %Ynew
            ProjmY=Yold*(mybetamat(:,1:j-1)/(mybetamat(:,1:j-1)'* mybetamat(:,1:j-1))* mybetamat(:,1:j-1)'); 
            ProjmY(isnan(ProjmY))=0;
            Y=Yold-ProjmY;    
       end
        
        while ( (iter<maxiteration) && (and ((diff_alpha > tol) , (diff_beta >tol) ))) %iterate if iteration is less than maxiteration or convergence not reached

            iter=iter+1;

            mycorr2=[];
            Tauy=median(Tauyvec);
            for jj=1:length(Tauxvec);
                 mycorr=[];
                 Taux=(Tauxvec(:,jj));
                  for ii=1:nfolds
                    which=foldid==ii;
                    [tildeA, tildeB, tilderho,~,~,~,~,~]=mynonsccaOptA(Xold(~which,:),Yold(~which,:)); 
                    if(iter==1)
                       myalphaold=tildeA(:,j);
                       mybetaold=tildeB(:,j);
                       tilderhoold=tilderho(j);
                    end
                    [myhatalpha, myhatbeta,myrho]= structuredsccaOptA(X(~which,:),Y(~which,:),mybetaold, myalphaold, tilderhoold,Taux,Tauy,edgesX,edgesY,weightsX,weightsY,mygamma,myeta,method);    
                    
                    %apply estimated canonical vectors on testing folds
                     myUtest=X(which,:)*myhatalpha;
                     myVtest=Y(which,:)*myhatbeta;

                    %apply estimated canonical vectors on training folds
                     myUtrain=X(~which,:)*myhatalpha;
                     myVtrain=Y(~which,:)*myhatbeta;

                     mycorrTrain=abs(corr(myUtrain,myVtrain));
                     mycorrTrain(isnan(mycorrTrain))=0;
                     mycorrTest=abs(corr(myUtest,myVtest));
                     mycorrTest(isnan(mycorrTest))=0;

                     mycorr=[mycorr;[Taux' ii mycorrTrain  mycorrTest ]];

                  end
                 isnanTrain=isnan(mycorr(:,3));
                 isnanTest=isnan(mycorr(:,4));    
             mycorr2=[mycorr2;[Taux'  mean(mycorr(~isnanTrain,3)) mean(mycorr(~isnanTest,4))]];
            end

            %choose lambda yielding minimum mean difference
            mydiffmean=abs(mycorr2(:,2)- mycorr2(:,3));
            if(all(mydiffmean)==0)
               row=[];
            else
             row= find(mydiffmean==nanmin(mydiffmean(mydiffmean~=0)),1,'last');
            end
             optTaux=Tauxvec(:,row); 
             if(isempty(row))
                optTaux=Tauxvec(:,1);
             end             
            
            mycorr2=[];
            Taux=optTaux;
            for jj=1:length(Tauyvec);
                mycorr=[];
                 Tauy=(Tauyvec(:,jj));
                   for ii=1:nfolds
                       which=foldid==ii;
                      [tildeA, tildeB, tilderho,~,~,~,~,~]=mynonsccaOptA(Xold(~which,:),Yold(~which,:));              
                       if(iter==1)
                           myalphaold=tildeA(:,j);
                           mybetaold=tildeB(:,j);
                           tilderhoold=tilderho(j);
                       end
                     [myhatalpha, myhatbeta,myrho]= structuredsccaOptA(X(~which,:),Y(~which,:),mybetaold, myalphaold, tilderhoold,Taux,Tauy,edgesX,edgesY,weightsX,weightsY,mygamma,myeta,method);

                    %apply estimated canonical vectors on testing folds
                     myUtest=Xold(which,:)*myhatalpha;
                     myVtest=Yold(which,:)*myhatbeta;

                    %apply estimated canonical vectors on training folds
                     myUtrain=Xold(~which,:)*myhatalpha;
                     myVtrain=Yold(~which,:)*myhatbeta;

                     mycorrTrain=abs(corr(myUtrain,myVtrain));
                     mycorrTrain(isnan(mycorrTrain))=0;
                     mycorrTest=abs(corr(myUtest,myVtest));
                     mycorrTest(isnan(mycorrTest))=0;

                     mycorr=[mycorr;[Tauy' ii mycorrTrain  mycorrTest ]];
                   end
                  isnanTrain=isnan(mycorr(:,3));
                 isnanTest=isnan(mycorr(:,4));    
             mycorr2=[mycorr2;[Tauy'  mean(mycorr(~isnanTrain,3)) mean(mycorr(~isnanTest,4))]];
            end

            %choose lambda yielding minimum mean difference
            mydiffmean=abs(mycorr2(:,2)-mycorr2(:,3));
            if (all(mydiffmean)==0)
              row=[];
            else
             row= find(mydiffmean==nanmin(mydiffmean(mydiffmean~=0)),1,'last');
            end
            optTauy=Tauyvec(row); 
             if(isempty(row))
                optTauy=Tauyvec(:,1);
             end   
             %apply chosen lambda on training set to update alpha and beta
            [tildeA, tildeB, tilderho,~,~,~,~,~]=mynonsccaOptA(Xold,Yold);      

            myalphaold=tildeA(:,j);
            mybetaold=tildeB(:,j);
            tilderhoold=tilderho(j);
                if(iter>1);
                    myalphaold=myalpha;
                    mybetaold=mybeta;
                    tilderhoold=mytilderho;
                end
            [myalpha, mybeta,mytilderho]=structuredsccaOptA(X,Y,mybetaold, myalphaold,...
                tilderhoold,optTaux,optTauy,edgesX,edgesY,weightsX,weightsY,mygamma,myeta,method);

             myoptTau=[myoptTau;[iter j optTaux optTauy]];

             if( or( min(sum(abs(myalpha)))==0, min(sum(abs(mybeta)))==0) )
                 myalpha=myalphaold;
                 mybeta=mybetaold;
                 break
             end  
             if (abs(tilderhoold-mytilderho) <= 10^-3)
                 myalpha=myalphaold;
                 mybeta=mybetaold;
                 break
             end
                 
            %convergence criterion and update
             diff_alpha=norm(abs(myalpha)-abs(myalphaold),Inf);
             diff_beta=norm(abs(mybeta)-abs(mybetaold),Inf);
             myalphaold=myalpha;
             mybetaold=mybeta;
        end
        
myest_corr(j,1)=abs(corr(Xold*myalpha, Yold*mybeta));
myalphamat(:,j)=myalpha;
mybetamat(:,j)=mybeta;
               
end

%display some results
NonZeroAlpha=sum(myalphamat~=0);
NonZeroBeta=sum(mybetamat~=0);
dNZalpha=['Num non-zeros alpha: ', num2str(NonZeroAlpha)];
disp(dNZalpha);
dNZbeta=['Num non-zeros beta: ', num2str(NonZeroBeta)];
disp(dNZbeta);
dcorr=['Corr(X\alpha,Y\beta): ', num2str(myest_corr')];
disp(dcorr);
d=['Structured Sparse CCA Method used is: ', method, ' Option A'];
disp(d);

%output results    
cvout.optTauX=myoptTau(:,1:3);%[iteration ncancorr optTayX];
cvout.optTauY=myoptTau(:,[1:2,4]);%[iteration ncancorr optTayY];
cvout.hatalpha=myalphamat;
cvout.hatbeta=mybetamat;
cvout.maxcorr=myest_corr;
cvout.method=method;
cvout.tunerange=[Tauxvec Tauyvec];
%%
%Author: Sandra E. Safo
%April 15, 2017
%structured sparse CCA examples for Option A
%load example data
clearvars; clc;
rng('default')

% This function depends on the CVX package. Please download at
%http://cvxr.com/cvx/download/

%addpath(....);

disp('  ');
disp('--------------------------------------');
disp('Example One Data'); 
disp('n = 80, p = q = 500');
disp('One CCA Vector');

disp('--------------------------------------');
disp('loading example one data ...');


load exampleonedata.mat;

%call cvstructuredscca to estimate optimal tuning parameter using training
%data;
%cvout=cvmultistructuredsccaOptA(X,Y,ncancorr,edgesX,edgesY,weightsX,weightsY,method,mygamma,myeta,nfolds,ngrid);

%required parameters
method='Fused' ; %'Fused' or 'Grouped';
ncancorr=1;
%will use default settings
tic;
cvoutGA=cvmultistructuredsccaOptA(X,Y,ncancorr,edgesX,edgesY,weightsX,weightsY,method);
ttGA=toc;
%Alphaerr = norm(cvout.hatalpha * cvout.hatalpha'  - mytruealpha * mytruealpha', 'fro');
%Betaerr = norm(cvout.hatbeta * cvout.hatalpha'  - mytruealpha * mytruealpha', 'fro');

%apply optimal tuning parameter on testing data
TauX=cvoutGA.optTauX;
TauY=cvoutGA.optTauY;
[maxcorr,hatalpha,hatbeta]= multistructuredsccaOptA(Xtest,Ytest,ncancorr,TauX,TauY,edgesX,edgesY,weightsX,weightsY,method);



%plot first CCA vectors with training data
figure()
scoresX=X*cvout.hatalpha;
scoresY=Y*cvout.hatbeta;
h1 = scatter(scoresX,scoresY,50,'r','filled');
set(h1,'LineWidth',2)
set(gca,'Fontsize',14);
title('Scores of X vs Y, Training Data');
box on;
set(gca,'LineWidth',2)

%plot first CCA vectors with testing data
figure()
scoresX=Xtest*hatalpha;
scoresY=Ytest*hatbeta;
h1 = scatter(scoresX,scoresY,50,'filled');
set(h1,'LineWidth',2)
set(gca,'Fontsize',14);
title('Scores of X vs Y, Testing Data');
box on;
set(gca,'LineWidth',2)

%%
%example two data
disp('  ');
disp('--------------------------------------');
disp('Example Two Data'); 
disp('n = 80, p = q = 500');
disp('Two CCA Vectors');

disp('--------------------------------------');
disp('loading example two data ...');


load exampletwodata.mat;

%call cvstructuredscca to estimate optimal tuning parameter using training
%data;
%cvout=cvmultistructuredsccaOptA(X,Y,ncancorr,edgesX,edgesY,weightsX,weightsY,method,mygamma,myeta,nfolds,ngrid);

%required parameters
method='Fused' ; %'Fused' or 'Grouped';
ncancorr=2;
%will use default settings
cvout=cvmultistructuredsccaOptA(X,Y,ncancorr,edgesX,edgesY,weightsX,weightsY,method);
%Alphaerr = norm(cvout.hatalpha * cvout.hatalpha'  - mytruealpha * mytruealpha', 'fro');
%Betaerr = norm(cvout.hatbeta * cvout.hatalpha'  - mytruealpha * mytruealpha', 'fro');

%apply optimal tuning parameter on testing data
TauX=cvout.optTauX;
TauY=cvout.optTauY;
[maxcorr,hatalpha,hatbeta]= multistructuredsccaOptA(Xtest,Ytest,ncancorr,TauX,TauY,edgesX,edgesY,weightsX,weightsY,method);



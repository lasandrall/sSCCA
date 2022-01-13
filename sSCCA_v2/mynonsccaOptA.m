function[tildeA, tildeB, tilderho,S12mybetaold,S21myalphaold,Ux,Uy,Sigmaxyr]=mynonsccaOptA(X,Y)
%function to obtain nonsparse solution to cca problem

%%%%Inputs
%X:         n x p dataset with network information, where rows are samples, and
%           columns are variables
%Y:         n x q dataset without network information, where rows are samples,
%           and columns are variables

[n,p]=size(X);
[n,q]=size(Y);
     
    X=X';Y=Y';
    [Ux, Dx, Vx] = svd(X, 'econ');
    Rx = Dx*Vx';
    cRx=Rx-repmat(mean(Rx,2),1,n);

    [Uy, Dy, Vy] = svd(Y, 'econ');
    Ry = Dy*Vy';
    cRy=Ry-repmat(mean(Ry,2),1,n);
    Sigmaxyr=(cRx*cRy')/(n-1);

    
    lambdax=sqrt(log(p)/n);
    lambday=sqrt(log(q)/n);

    Sigmaxxr=(cRx*cRx')/(n-1);
    Sigmayyr=(cRy*cRy')/(n-1);
    [tildeAr, Dr, tildeBr]=svd( sqrtm(inv(Sigmaxxr +lambdax*eye(n,n) ))*Sigmaxyr*sqrtm(inv(Sigmayyr +lambday*eye(n,n) )));
    %Sigmaxy=Ux*Sigmaxyr*Uy';

    tilderho=diag(Dr);

    tildeA=Ux*sqrtm(inv(Sigmaxxr+ lambdax*eye(n,n)))*tildeAr;
    tildeA=tildeA./repmat(sum((tildeA.*tildeA),1).^(0.5),p,1);

    tildeB=Uy*sqrtm(inv(Sigmayyr+lambday*eye(n,n)))*tildeBr;
    tildeB=tildeB./repmat(sum((tildeB.*tildeB),1).^(0.5),q,1);

    S12mybetaold=Ux*Sigmaxyr*Uy'*tildeB;
    S21myalphaold=Uy*Sigmaxyr'*Ux'*tildeA;



%Sandra E. Safo
%All rights reserved
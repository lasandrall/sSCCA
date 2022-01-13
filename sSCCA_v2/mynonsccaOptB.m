function[tildeA, tildeB, tilderho,SigxxinvSigxy,SigyyinvSigyx]=mynonsccaOptB(X,Y)
%function to obtain nonsparse solution to cca problem

%%%%Inputs
%X:         n x p dataset with network information, where rows are samples, and
%           columns are variables
%Y:         n x q dataset without network information, where rows are samples,
%           and columns are variables
%%%%Outputs
%tildeA:    nonsparse CCA solution for X
%tildeB:    nonsparse CCA solution for Y
%tilderho:  Canonical correlaiton coefficient
%SigxxinvSigxy and SigyyinvSigyx: matrices for the constraint functions
    [n,p]=size(X);
    [n,q]=size(Y);
   
if(or(p>n,q>n)) %high dimensional setting    
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
    Sigmaxy=Ux*Sigmaxyr*Uy';

    tilderho=diag(Dr);

    tildeA=Ux*sqrtm(inv(Sigmaxxr+lambdax*eye(n,n)))*tildeAr;
    tildeA=tildeA./repmat(sum((tildeA.*tildeA),1).^(0.5),p,1);

    tildeB=Uy*sqrtm(inv(Sigmayyr+lambday*eye(n,n)))*tildeBr;
    tildeB=tildeB./repmat(sum((tildeB.*tildeB),1).^(0.5),q,1);

   
    SigxxinvSigxy=Ux*inv(Sigmaxxr +lambdax*eye(n,n))*Sigmaxyr*Uy';
    SigyyinvSigyx=Uy*inv(Sigmayyr +lambday*eye(n,n))*Sigmaxyr'*Ux';
    
else
cx=X-repmat(mean(X,1),n,1);

cy=Y-repmat(mean(Y,1),n,1);
Sigmaxy=(cx'*cy)/(n-1);

Sigmaxx=(cx'*cx)/(n-1);
Sigmayy=(cy'*cy)/(n-1);
[tildeA, D, tildeB]=svd( sqrtm(inv(Sigmaxx ))*Sigmaxy*sqrtm(inv(Sigmayy)));

tildeA=sqrtm(inv(Sigmaxx))*real(tildeA);
tildeB=sqrtm(inv(Sigmayy))*real(tildeB);
tilderho=real(diag(D));

SigxxinvSigxy=Sigmaxx\Sigmaxy;
SigyyinvSigyx=Sigmayy\Sigmaxy';
end

%Sandra E. Safo
%All rights reserved
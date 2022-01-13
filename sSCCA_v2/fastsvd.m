function [U,Sigmar]=fastsvd(mydata);

%mydata is a n by p matrix
    [n,p]=size(mydata);
       
    
    mydata=mydata';
    [U, D, V] = svd(mydata, 'econ');
    R = D*V';
    cR=R-repmat(mean(R,2),1,n);


    Sigmar=(cR*cR')/(n-1);
   

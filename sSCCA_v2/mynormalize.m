%my normalized data
function X=mynormalize(mydata);
[N,p]=size(mydata);
mycenter=mydata-repmat(mean(mydata,1),N,1);
mystddata=mycenter./repmat(std(mydata,0),N,1);

X=mystddata;
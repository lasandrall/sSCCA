function [myalpha, mybeta,myrho]= structuredsccaOptA(X,Y,mybetaold, myalphaold, tilderhoold,Taux,Tauy,edgesX,edgesY,weightsX,weightsY,mygamma,eta,method)
%this calls the CVX algorithm to estimate alpha and beta

[~, ~, ~,~,~,Ux,Uy,Sigmaxyr]=mynonsccaOptA(X,Y);

p=size(X,2);
[n,q]=size(Y);

[~,idx]=sort(edgesX(:,1));
edges2X=edgesX(idx,:);

[~,idy]=sort(edgesY(:,1));
edges2Y=edgesY(idy,:);

      
        N=size(edgesX,1);
        cvx_begin quiet
                variable alphai(p)
                alphai_w = alphai./ weightsX.^ (1/mygamma) ;
                if(strcmp(method,'Grouped'))
                sumedges=sum(norms(alphai_w(edges2X),mygamma,2));
                elseif(strcmp(method,'Fused'))
                sumedges=sum(abs(alphai(edges2X(:,1))./weightsX(edges2X(:,1))-alphai(edges2X(:,2))./weightsX(edges2X(:,2)))); 
                end
                minimize((1-eta)*sumedges + eta*norm(alphai(setdiff(1:p,unique(edgesX))),1))
                [U,Sigmar]=fastsvd(X);
                Utralphai=U'*alphai;
                subject to
                Sxymybetaold=Ux*Sigmaxyr*Uy'*mybetaold;
                norm(Sxymybetaold(:,1)-tilderhoold(1)*(U*Sigmar*Utralphai + sqrt(log(p)/n)*alphai),Inf)<=Taux;
          cvx_end

        N=size(edgesY,1);
        cvx_begin quiet
                variable betai(q)
                betai_w = betai./ weightsY.^ (1/mygamma) ;
                if(strcmp(method,'Grouped'))
                sumedges=sum(norms(betai_w(edges2Y),mygamma,2));
                elseif(strcmp(method,'Fused'))
                 sumedges=sum(abs(betai(edges2Y(:,1))./weightsY(edges2Y(:,1))-betai(edges2Y(:,2))./weightsY(edges2Y(:,2)))); 
                end
                minimize((1-eta)*sumedges + eta*norm(betai(setdiff(1:q,unique(edgesY))),1));
                [U,Sigmar]=fastsvd(Y);
                Utrbetai=U'*betai;
                subject to
                Syxmyalphaold=Uy*Sigmaxyr'*Ux'*myalphaold;
                norm(Syxmyalphaold(:,1)-tilderhoold(1)*(U*Sigmar*Utrbetai + sqrt(log(q)/n)*betai),Inf)<=Tauy;
        cvx_end
        %normalize alpha
        alphai(abs(alphai)<=10^(-5))=0;
        if(sum(abs(alphai))~=0)
            myalpha=alphai./norm(alphai,2);
        else
            myalpha=alphai;
        end
         

        %normalize beta
        betai(abs(betai)<=10^(-5))=0;
        if(sum(abs(betai))~=0)
         mybeta=betai./norm(betai,2);
        else
            mybeta=betai;
        end

        if(or(sum(abs(myalpha))==0, sum(abs(mybeta))==0))
         myrho=0;
        else
         myrho=abs(corr(X*myalpha, Y*mybeta));
        end




        
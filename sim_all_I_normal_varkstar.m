clear all;
clc
diary('sim_all_I_normal_varkstar.txt');
%diary('sim_all_t3.txt');
%diary('sim_all_chi3.txt');
pp=parpool('local',25);
tic
N1=1000;
N2=40;
res_table=cell(2,2);
n_ind=[100,500,1000,2000];



for alphafor=1:2
    
    alpha=0.2*alphafor;

    for cfor=1:2
        c=0.2*cfor;
        
        for t=1:4
            
            n=n_ind(t);
            p=round(c*n);
            k=round(alpha*n);
            k0=round(k/2);
            para_theta0=1;
            
            X=unifrnd(1,5,n,k); 

            X0=X(:,1:k0);
            theta0=zeros(p,1);
            parfor i=1:p
                theta0(i,1)=(-1/2)^(i-1);
            end
            theta0=para_theta0*theta0;
            Theta0=ones(k0,1)*theta0';
            Theta=[Theta0;zeros(k-k0,p)];
            P=X/(X'*X)*X';
            Q=eye(n)-P;
            A=zeros(n,k);
            parfor j=1:k
                Xj=X;
                Xj(:,j)=[];
                Qj=eye(n)-Xj/(Xj'*Xj)*Xj';
                A(:,j)=Qj*X(:,j)/sqrt(X(:,j)'*Qj*X(:,j));
            end
            maxZ=zeros(1,N1);
            %%%%%%%%%%%%% Case I %%%%%%%%%%%%%%%%%%%%%
            parfor mm=1:N1
                Z1=randn(n,p);
                maxZ(1,mm)=max(diag(A'*Z1/(Z1'*Q*Z1/n)*Z1'*A/n));
            end
       
            

            maxZ0=max(maxZ);
            maxZ5=prctile(maxZ,95);
            maxZ10=prctile(maxZ,90);
            maxZ50=median(maxZ);
            
            res_cr=zeros(5,3,N2);
            res_sz=zeros(5,3,N2);
            
            spmd
                for i=1:N2
                    E=randn(n,p); %Normal distribution
                    %E=trnd(3,n,p)/sqrt(3); %t3 distribution
                    %E=(chi2rnd(3,n,p)-3)/sqrt(6); %chi3 distribution
                    Y=X*Theta+E;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    K=diag(A'*Y/(E'*Q*E/n)*Y'*A/n);
                    AIC=log(1+K)-2*c;
                    BIC=log(1+K)-log(n)*c;
                    Cp=(1-alpha)*K-2*c;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    truemodel=[1:k0]';
                    
                    AICj_star=find(AIC>0);
                    AIC_sz=length(AICj_star);
                    if isequal(truemodel,AICj_star)
                        res_cr(1,2,i)=1;
                        res_sz(1,2,i)=AIC_sz;
                    elseif all(ismember(truemodel,AICj_star))
                        res_cr(1,3,i)=1;
                        res_sz(1,3,i)=AIC_sz;
                    else
                        res_cr(1,1,i)=1;
                        res_sz(1,1,i)=AIC_sz;
                    end
                    
                    
                    BICj_star=find(BIC>0);
                    BIC_sz=length(BICj_star);
                    if isequal(truemodel,BICj_star)
                        res_cr(2,2,i)=1;
                        res_sz(2,2,i)=BIC_sz;
                    elseif all(ismember(truemodel,BICj_star))
                        res_cr(2,3,i)=1;
                        res_sz(2,3,i)=BIC_sz;
                    else
                        res_cr(2,1,i)=1;
                        res_sz(2,1,i)=BIC_sz;
                    end
                    
                    
                    Cpj_star=find(Cp>0);
                    Cp_sz=length(Cpj_star);
                    if isequal(truemodel,Cpj_star)
                        res_cr(3,2,i)=1;
                        res_sz(3,2,i)=Cp_sz;
                    elseif all(ismember(truemodel,Cpj_star))
                        res_cr(3,3,i)=1;
                        res_sz(3,3,i)=Cp_sz;
                    else
                        res_cr(3,1,i)=1;
                        res_sz(3,1,i)=Cp_sz;
                    end
                    
                    
                    KOO0_star=find(K>maxZ0);
                    KOO0_sz=length(KOO0_star);
                    if isequal(truemodel,KOO0_star)
                        res_cr(4,2,i)=1;
                        res_sz(4,2,i)=KOO0_sz;
                    elseif all(ismember(truemodel,KOO0_star))
                        res_cr(4,3,i)=1;
                        res_sz(4,3,i)=KOO0_sz;
                    else
                        res_cr(4,1,i)=1;
                        res_sz(4,1,i)=KOO0_sz;
                    end
                    
                    KOO5_star=find(K>maxZ5);
                    KOO5_sz=length(KOO5_star);
                    if isequal(truemodel,KOO5_star)
                        res_cr(5,2,i)=1;
                        res_sz(5,2,i)=KOO5_sz;
                    elseif all(ismember(truemodel,KOO5_star))
                        res_cr(5,3,i)=1;
                        res_sz(5,3,i)=KOO5_sz;
                    else
                        res_cr(5,1,i)=1;
                        res_sz(5,1,i)=KOO5_sz;
                    end
                    
%                     KOO10_star=find(K>maxZ10);
%                     KOO10_sz=length(KOO10_star);
%                     if isequal(truemodel,KOO10_star)
%                         res_cr(6,2,i)=1;
%                         res_sz(6,2,i)=KOO10_sz;
%                     elseif all(ismember(truemodel,KOO10_star))
%                         res_cr(6,3,i)=1;
%                         res_sz(6,3,i)=KOO10_sz;
%                     else
%                         res_cr(6,1,i)=1;
%                         res_sz(6,1,i)=KOO10_sz;
%                     end
%                     
%                     KOO50_star=find(K>maxZ50);
%                     KOO50_sz=length(KOO50_star);
%                     if isequal(truemodel,KOO50_star)
%                         res_cr(7,2,i)=1;
%                         res_sz(7,2,i)=KOO50_sz;
%                     elseif all(ismember(truemodel,KOO50_star))
%                         res_cr(7,3,i)=1;
%                         res_sz(7,3,i)=KOO50_sz;
%                     else
%                         res_cr(7,1,i)=1;
%                         res_sz(7,1,i)=KOO50_sz;
%                     end                    
                end
                sum_cr=sum(res_cr,3);
                sum_sz=sum(res_sz,3);
            end
            
            parameter=['The value of n c alpha are ' num2str([n c alpha])];
            disp(parameter);
            
            sumres_cr=sum(cat(3, sum_cr{:}),3)
            sumres_sz=sum(cat(3, sum_sz{:}),3)./sumres_cr
           
            sumres=[sumres_cr,sumres_sz(:,3)-k0];

            res_table{alphafor,cfor}(:,:,t)=sumres;
            
            toc
        end
        save sim_all_I_normal_varkstar.mat res_table
        %save sim_all_t3.mat plotres_cr plotres_sz
        %save sim_all_chi3.mat plotres_cr plotres_sz
    end
end
delete(pp);


clear input
Rname={'Under';'True';'Over';'Average'};
input.data=table([res_table{1,1}(:,:,3)',res_table{1,1}(:,:,4)'],'RowNames',Rname);
input.dataFormat = {'%i'};
input.tableColumnAlignment = 'c';
input.tableBorders = 1;
input.makeCompleteLatexDocument = 1;
latex = latexTable(input);

clear input
input.data=[res_table{1,1}(:,:,1)',res_table{1,1}(:,:,2)';...
    res_table{1,1}(:,:,3)',res_table{1,1}(:,:,4)';
    res_table{1,2}(:,:,1)',res_table{1,2}(:,:,2)';...
    res_table{1,2}(:,:,3)',res_table{1,2}(:,:,4)';
    res_table{2,1}(:,:,1)',res_table{2,1}(:,:,2)';...
    res_table{2,1}(:,:,3)',res_table{2,1}(:,:,4)';
    res_table{2,2}(:,:,1)',res_table{2,2}(:,:,2)';...
    res_table{2,2}(:,:,3)',res_table{2,2}(:,:,4)'];
input.tableRowLabels={'Under';'True';'Over';'Average';...
    'Under';'True';'Over';'Average';...
    'Under';'True';'Over';'Average';...
    'Under';'True';'Over';'Average';...
    'Under';'True';'Over';'Average';...
    'Under';'True';'Over';'Average';...
    'Under';'True';'Over';'Average';...
    'Under';'True';'Over';'Average'}
input.dataFormat = {'%.2f'};
input.dataNanString = '--';
input.tableColumnAlignment = 'c';
input.tableBorders = 1;
input.makeCompleteLatexDocument = 1;
latex = latexTable(input);
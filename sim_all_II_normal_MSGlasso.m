clear;
clc
diary('sim_all_II_normal_MSGlasso.txt');
%diary('sim_all_t3.txt');
%diary('sim_all_chi3.txt');
%pp=parpool('local',25);
tic
N1=100;
N2=40;
plotres_cr=cell(2,2);
plotres_sz=cell(2,2);
n_ind=[100,500,1000,2000];
k0=5;


for alphafor=1:2

    alpha=0.2*alphafor;

    for cfor=1:2
        c=0.2*cfor;

        for t=1:3

            n=n_ind(t);
            p=round(c*n);
            k=round(alpha*n);
            para_theta0=1;

            X=unifrnd(1,5,n,k);

            X0=X(:,1:k0);
            theta0=zeros(p,1);
            for i=1:p
                theta0(i,1)=(-1/2)^(i-1);
            end
            theta0=para_theta0*theta0;
            Theta0=ones(k0,1)*theta0';
            Theta=[Theta0;zeros(k-k0,p)];
            P=X/(X'*X)*X';
            Q=eye(n)-P;
            A=zeros(n,k);
            for j=1:k
                Xj=X;
                Xj(:,j)=[];
                Qj=eye(n)-Xj/(Xj'*Xj)*Xj';
                A(:,j)=Qj*X(:,j)/sqrt(X(:,j)'*Qj*X(:,j));
            end
            maxZ=zeros(1,N1);
            %%%%%%%%%%%%% Case I %%%%%%%%%%%%%%%%%%%%%
            for mm=1:N1
                Z1=randn(n,p);
                maxZ(1,mm)=max(diag(A'*Z1/(Z1'*Q*Z1/n)*Z1'*A/n));
            end



            maxZ0=max(maxZ);
            maxZ5=prctile(maxZ,95);
            maxZ10=prctile(maxZ,90);
            maxZ50=median(maxZ);

            res_cr=zeros(5,3,N2);
            res_sz=zeros(5,3,N2);

            L=diag([1:p]);
            [QQ,RR] = qr(unifrnd(-5,5,p,p));

            Sigma=QQ*L*QQ';
            sqrtSigma=QQ*sqrt(L)*QQ';
            E=randn(n,p); %Normal distribution
            %E=trnd(3,n,p)/sqrt(3); %t3 distribution
            %E=(chi2rnd(3,n,p)-3)/sqrt(6); %chi3 distribution
            Y=X*Theta*sqrtSigma+E*sqrtSigma;
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            save XY4lam1G.mat X Y
% 
%             Rpath = '/usr/local/bin/';
%             lam1G = '/Users/john/Library/CloudStorage/Dropbox/huj/Fujikoshi/model_selection/my_version/Sinica/revision/simulation/simall-20240130/lam1G.R';
%             MSGlasso= '/Users/john/Library/CloudStorage/Dropbox/huj/Fujikoshi/model_selection/my_version/Sinica/revision/simulation/simall-20240130/sim_all_MSGlasso.R';
%            RunRcode_mac(lam1G, Rpath);

            Rpath = 'C:/Program Files/R/R-4.3.2/bin';
            lam1G = 'C:/Users/admin/OneDrive/simall-20240130/lam1G.R';
            MSGlasso= 'C:/Users/admin/OneDrive/simall-20240130/sim_all_MSGlasso.R';
            RunRcode(lam1G, Rpath);

            truemodel=[1:k0]';

            for i=1:N1
                E=randn(n,p)*sqrtSigma; %Normal distribution
                %E=trnd(3,n,p)/sqrt(3); %t3 distribution
                %E=(chi2rnd(3,n,p)-3)/sqrt(6); %chi3 distribution
                Y=X*Theta*sqrtSigma+E;
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                K=diag(A'*Y/(E'*Q*E/n)*Y'*A/n);
                %%%%%%%%%%%%%%%%%%%%%%%%%%

                save modelXY.mat X Y
                RunRcode(MSGlasso, Rpath);
                Thetatem = load('hatTheta.mat');
                delete('hatTheta.mat');
                hatTheta=Thetatem.hattheta;

                MSG = max(hatTheta,[],2);

                MSGj_star=find(MSG>0.00000000000000001);
                MSG_sz=length(MSGj_star);
                if isequal(truemodel,MSGj_star)
                    res_cr(1,2,i)=1;
                    res_sz(1,2,i)=MSG_sz;
                elseif all(ismember(truemodel,MSGj_star))
                    res_cr(1,3,i)=1;
                    res_sz(1,3,i)=MSG_sz;
                else
                    res_cr(1,1,i)=1;
                    res_sz(1,1,i)=MSG_sz;
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


            end
            sum_cr=sum(res_cr,3);
            sum_sz=sum(res_sz,3);


            parameter=['The value of n c alpha are ' num2str([n c alpha])];
            disp(parameter);

            sumres_cr=sum(res_cr,3)
            sumres_sz=sum(res_sz,3)./sumres_cr

            plotres_cr{alphafor,cfor}(:,:,t)=sumres_cr;
            plotres_sz{alphafor,cfor}(:,:,t)=sumres_sz;
            toc
        end

        save plotres_II_normal.mat plotres_cr plotres_sz
        %save sim_all_t3.mat plotres_cr plotres_sz
        %save sim_all_chi3.mat plotres_cr plotres_sz
    end
end




% clear input
% Rname={'Under';'True';'Over';'Average'};
% input.data=table([res_table{1,1}(:,:,3)',res_table{1,1}(:,:,4)'],'RowNames',Rname);
% input.dataFormat = {'%i'};
% input.tableColumnAlignment = 'c';
% input.tableBorders = 1;
% input.makeCompleteLatexDocument = 1;
% latex = latexTable(input);

% clear input
% input.data=[res_table{1,1}(:,:,1)',res_table{1,1}(:,:,2)';...
%     res_table{1,1}(:,:,3)',res_table{1,1}(:,:,4)';
%     res_table{1,2}(:,:,1)',res_table{1,2}(:,:,2)';...
%     res_table{1,2}(:,:,3)',res_table{1,2}(:,:,4)';
%     res_table{2,1}(:,:,1)',res_table{2,1}(:,:,2)';...
%     res_table{2,1}(:,:,3)',res_table{2,1}(:,:,4)';
%     res_table{2,2}(:,:,1)',res_table{2,2}(:,:,2)';...
%     res_table{2,2}(:,:,3)',res_table{2,2}(:,:,4)'];
% input.tableRowLabels={'Under';'True';'Over';'Average';...
%     'Under';'True';'Over';'Average';...
%     'Under';'True';'Over';'Average';...
%     'Under';'True';'Over';'Average';...
%     'Under';'True';'Over';'Average';...
%     'Under';'True';'Over';'Average';...
%     'Under';'True';'Over';'Average';...
%     'Under';'True';'Over';'Average'}
% input.dataFormat = {'%.2f'};
% input.dataNanString = '--';
% input.tableColumnAlignment = 'c';
% input.tableBorders = 1;
% input.makeCompleteLatexDocument = 1;
% latex = latexTable(input);
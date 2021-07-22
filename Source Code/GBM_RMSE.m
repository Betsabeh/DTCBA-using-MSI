function [AUPR,RMSE,CI,RM2]=GBM_RMSE(X_KD,X_S,KD_Bind,X_test_KD,X_test_S,L_KD_test)
% This function run GBM to combine and predict binding affinity values
%--------------------------------------------------------------------------
% [AUPR,RMSE,CI,RM2]=GBM_RMSE(X_KD,X_S,KD_Bind,X_test_KD,X_test_S,L_KD_test)
% Inputs:
% X_KD:       matrix of KNN binding affinity values for tarin data (N_train*25*KNN)
% X_S:        matrix of KNN similarity values for train data (N_train*25*KNN)
% KD_Bind:    Binding affinity of train data
% X_test_KD:  matrix of KNN binding affinity values for test data (N_test*25*KNN)
% X_test_S:   matrix of KNN similarity values for test data (N_test*25*KNN)
% L_KD_test:  Binding affinity of test data
%
% Outputs:
% AUPR,RMSE,CI,RM2
%--------------------------------------------------------------------------
global Bag_rate
global max_iter
global shrink_rate
global max_leaf
global KNN

num_Bind=numel(KD_Bind);
num_f=size(X_KD,2);
num_test=numel(L_KD_test);

avg=mean(KD_Bind);
F_Bind=ones(num_Bind,1).*avg;
F_test=ones(num_test,1)*avg;
num_tr=ceil(Bag_rate*num_Bind);

for i=1:max_iter
    r=randi(num_f,1);
    feat_index{i}=randperm(num_f,r);
end
it=1;
J(1,1)=0;
l=1;
%--------------------------------------------------------------------------
while (it<=max_iter)
    % Gradient J1 RMSE
    temp=(KD_Bind-F_Bind);
    J(it)=sqrt(mean(temp.^2));
    Gradient_Bind=(temp./J(it));%(1/num_Bind)*
    %% Model
    Fix_S=randperm(num_Bind,num_tr);%ind;%
    X_Bind=[X_KD(Fix_S,feat_index{it},:),X_S(Fix_S,feat_index{it},:)];
    X_Bind=reshape(X_Bind,num_tr,2*numel(feat_index{it})*KNN);
    
    mdl=fitrtree(X_Bind,Gradient_Bind(Fix_S),'MaxNumSplits',max_leaf);
    % save(strcat('mdl_',num2str(it)),'mdl');
    X_Bind=[X_KD(:,feat_index{it},:),X_S(:,feat_index{it},:)];
    X_Bind=reshape(X_Bind,num_Bind,2*numel(feat_index{it})*KNN);
    temp_Bind=predict(mdl,X_Bind);
    [Learning_rate,~]=fminunc(@(Learning_rate)learning_rate(KD_Bind,F_Bind,temp_Bind,Learning_rate),0);%,options);
    L_rate=Learning_rate*shrink_rate;
    F_Bind=F_Bind+L_rate*temp_Bind;
    X_test=[X_test_KD(:,feat_index{it},:),X_test_S(:,feat_index{it},:)];
    X_test=reshape(X_test,num_test,2*numel(feat_index{it})*KNN);
    
    F_test=F_test+L_rate*predict(mdl,X_test);
    
    if mod(it,100)==0
        
      [AUPR(l),RMSE(l),CI(l),RM2(l)]=Validation(L_KD_test,F_test)
        
        
        l=l+1;
        
    end
    
   % fprintf('iteration=%d  MSE=%f\n',it,J(it))
    it=it+1;
end
% % % Test_AUC= calculate_auc (F_Bind,KD_Bind>-log10(th))
% % % Test_AUPR= calculate_aupr(F_Bind,KD_Bind>-log10(th))
% % % figure
% % % plot(J)
% % % J=J(100:100:max_iter);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function J=learning_rate(Y_Bind,F_Bind,temp_Bind,Learing_rate)
%this function try to find optimum learning rate for gradient boosting

F_Bind=F_Bind+Learing_rate.*temp_Bind;
temp=Y_Bind-F_Bind;
J=sqrt(mean((temp.^2)));
end
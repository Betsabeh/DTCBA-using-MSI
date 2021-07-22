function [AUPR,RMSE,CI,RM2]=fused_GB(KD_Train,ind_d_tr,ind_p_tr,KD_Test,ind_d_te,ind_p_te)
%This function creates pairwise similarities and combine them with GBM
%--------------------------------------------------------------------------
%[AUPR,RMSE,CI,RM2]=fused_GB(KD_Train,Interactions,ind_d_tr,ind_p_tr,KD_Test,ind_d_te,ind_p_te)
%
%  Inputs:
%  KD_Train:  Binding affinity of train data
%  KD_Test:   Binding affinity of test data
%  ind_d_tr:  train data, drug indexes related to the drug-drug similarities
%  ind_p_tr:  train data, protein indexes related to the target-target similarities
%  ind_d_te:  test data, drug indexes related to the drug-drug similarities
%  ind_p_te:  test data, protein indexes related to the target-target similarities
%
% Outputs:
% AUPR,RMSE,CI,RM2
%--------------------------------------------------------------------------
%% Paramters
N_train=numel(KD_Train);
N_test=numel(KD_Test);

%--------------------------------------------------------------------------
%% Test Set
if N_test<7000
    [KD_hat_te,S_te]=Find_KNN_features(KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,'Test');
else %% faseter for large set
   index=0:2000:N_test;
   index(end+1)=N_test;
    for i=1:numel(index)-1
    ind11=index(i)+1;
    ind12=index(i+1);
     [KD_hat_te(ind11:ind12,:,:),S_te(ind11:ind12,:,:)]=Find_KNN_features(KD_Train,ind_d_tr,ind_p_tr,ind_d_te(ind11:ind12),ind_p_te(ind11:ind12),'Test');
    end
end

%--------------------------------------------------------------------------
%% Train set
if N_train<25000
    [KD_hat_tr,S_tr]=Find_KNN_features(KD_Train,ind_d_tr,ind_p_tr,ind_d_tr,ind_p_tr,'Train');
else
    index=0:6000:N_train;
    index(end+1)=N_train;
    for i=1:numel(index)-1
    ind11=index(i)+1;
    ind12=index(i+1);
    [KD_hat_tr(ind11:ind12,:,:),S_tr(ind11:ind12,:,:)]=Find_KNN_features(KD_Train,ind_d_tr,ind_p_tr,ind_d_tr(ind11:ind12),ind_p_tr(ind11:ind12),'Train');
    end
end

%-------------------------------------------------------------------------
%% Optimization GBM
[AUPR,RMSE,CI,RM2]=GBM_RMSE(KD_hat_tr,S_tr,KD_Train,KD_hat_te,S_te,KD_Test);

end
%--------------------------------------------------------------------------


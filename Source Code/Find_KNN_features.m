function [Pred,S]=Find_KNN_features(KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,option)
% This function finds the KNN for pairwise similarity measures
%--------------------------------------------------------------------------
%[Pred,S]=Find_KNN_features(Interactions,KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,option)
% 
%Inputs:
%  KD_Train:  Binding affinity of train data
%  ind_d_tr:  train data, drug indexes related to the drug-drug similarities
%  ind_p_tr:  train data, protein indexes related to the target-target similarities
%  ind_d_te:  test data, drug indexes related to the drug-drug similarities
%  ind_p_te:  test data, protein indexes related to the target-target similarities
%  option:    determine one of three procedures (warm, NewDrug, NewProt)
%
%Outputs:
%   Pred:   K-nearest neighbor binding affinity values  (N*25*KNN)
%   S:      K-nearest neighbor similarities values      (N*25*KNN)
%--------------------------------------------------------------------------
global cv_setting
switch cv_setting
    case 'warm'
        [Pred,S]=TwoD_mode(KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,option);
    case  'NewDrug'
        [Pred,S]=Single_mode_D(KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,option);
    case  'NewProt'
        [Pred,S]=Single_mode_P(KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,option);
end

end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [Pred,S]=TwoD_mode(KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,option)
global S_D S_P KNN
if strcmp(option,'Train')==1
    NN_index=2:KNN+1;
else
    NN_index=1:KNN;
end
%--------
N1=numel(S_D);
N2=numel(S_P);
for i=1:numel(ind_d_te)
    h=1;
    index1=find(ind_d_tr==ind_d_te(i));
    index2=find(ind_p_tr==ind_p_te(i));
    ind=union(index1,index2);
    for l=1:N1
        temp_d=S_D{l};
        for j=1:N2
            temp_p=S_P{j};
            S_mat=temp_d(ind_d_te(i),ind_d_tr(ind)).*temp_p(ind_p_te(i),ind_p_tr(ind));
            [Val_p,index_p]=sort(S_mat,'descend');
            if numel(ind)<=KNN
                temp=NN_index;
                temp=temp(1):numel(ind);
                KD_NN=KD_Train(ind(index_p(temp)));
                total_Sim=Val_p(temp);
                Pred(i,h,1:numel(temp))=KD_NN';
                S(i,h,1:numel(temp))=total_Sim;
                
            else
                KD_NN=KD_Train(ind(index_p(NN_index)));
                total_Sim=Val_p(NN_index);
                Pred(i,h,1:KNN)=KD_NN';
                S(i,h,1:KNN)=total_Sim;
            end
            
            h=h+1;
        end
    end
end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [Pred,S]=Single_mode_D(KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,option)
%% for new drug situation
global S_D S_P KNN
if strcmp(option,'Train')==1
    NN_index=2:KNN+1;
else
    NN_index=1:KNN;
end
N1=4;%% we dont consider interation similarity
for i=1:numel(ind_d_te)
    h=1;
    ind=find(ind_p_tr==ind_p_te(i));
    for l=1:N1
        temp_d=S_D{l};
        % for j=1:N2
        %  temp_p=S_P{j};
        S_mat=temp_d(ind_d_te(i),ind_d_tr(ind));%.*temp_p(ind_p_te(i),ind_p_tr(ind));
        [Val_p,index_p]=sort(S_mat,'descend');
        if numel(ind)<=KNN
            temp=NN_index;
            temp=temp(1):numel(ind);
            KD_NN=KD_Train(ind(index_p(temp)));
            total_Sim=Val_p(temp);
            Pred(i,h,1:numel(temp))=KD_NN';
            S(i,h,1:numel(temp))=total_Sim;
            
        else
            KD_NN=KD_Train(ind(index_p(NN_index)));
            total_Sim=Val_p(NN_index);
            Pred(i,h,1:KNN)=KD_NN';
            S(i,h,1:KNN)=total_Sim;
        end
       
        h=h+1;
        %  end
    end
end

end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [Pred,S]=Single_mode_P(KD_Train,ind_d_tr,ind_p_tr,ind_d_te,ind_p_te,option)
%% for new target situation
global S_D S_P KNN
if strcmp(option,'Train')==1
    NN_index=2:KNN+1;
else
    NN_index=1:KNN;
end
N2=4;
for i=1:numel(ind_d_te)
    h=1;
    ind=find(ind_d_tr==ind_d_te(i));
    for j=1:N2
        temp_p=S_P{j};
        S_mat=temp_p(ind_p_te(i),ind_p_tr(ind));
        [Val_p,index_p]=sort(S_mat,'descend');
        if numel(ind)<=KNN
            temp=NN_index;
            temp=temp(1):numel(ind);
            KD_NN=KD_Train(ind(index_p(temp)));
            total_Sim=Val_p(temp);
            Pred(i,h,1:numel(temp))=KD_NN';
            S(i,h,1:numel(temp))=total_Sim;
            
        else
            KD_NN=KD_Train(ind(index_p(NN_index)));
            total_Sim=Val_p(NN_index);
            Pred(i,h,1:KNN)=KD_NN';
            S(i,h,1:KNN)=total_Sim;
        end
        
        h=h+1;
    end
end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


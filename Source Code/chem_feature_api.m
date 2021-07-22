function [Drug_feature_list,NOT_find]=chem_feature_api(Drug_pubchem_list)
%% this function find Drug finger print (Chemical Structure)
load('all_human_drug_id_api')
load('all_human_drug_feat_api')
n=numel(all_human_drug_id_api)+1;
[d_id, ~, ic] = unique(Drug_pubchem_list);
k=1;
NOT_find=[];
for i=1:numel(d_id)
    i
    index=find(ismember(all_human_drug_id_api,d_id{i})~=0);
    if isempty(index)==0
        Drug_feature_list(i,:)=all_human_drug_feat_api(index,:);
    else
        FP_code=pubchem_fingerprint_api(d_id{i});
        temp=int32(FP_code)-'0';
        all_human_drug_id_api{n}=d_id{i};
        all_human_drug_feat_api(n,:)=temp;
        n=n+1;
        Drug_feature_list(i,:)=temp;
        save('all_human_drug_id_api','all_human_drug_id_api');
        save('all_human_drug_feat_api','all_human_drug_feat_api');
        
    end
end

Drug_feature_list=Drug_feature_list(ic,:);
end

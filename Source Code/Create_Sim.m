function [ind_d,ind_p]=Create_Sim(Data)
%This function creates Similarity from different sources
%--------------------------------------------------------------------------
% [ind_d,ind_p]=Create_Sim(Data)
% Input:
% Data:   Dataset which contain drug and target properties like ( Pubchem
% id or Uniprot id or Gene name)
% Oputputs:
% ind_d: Drug index in the Drugs list
% ind_p: Protein index in the Prots list
%--------------------------------------------------------------------------

%Parameters
global Drugs Prots data_name
global S_D S_P
Drug_list=Data(:,2);
if strcmp(data_name,'Large_KIBA')==1  || strcmp(data_name,'Small_KIBA')==1  
    Prot_list=Data(:,5);  %KBIA for uniprot id =5
else
    Prot_list=Data(:,1); %%% Daivs and Metz 1 for gene ,
end
Smile_Seq=Data(:,4);
List_seq=Data(:,3);
[Drugs,index_d,~]=unique(Drug_list);
USmile=Smile_Seq(index_d,:);
[Prots,index_p,~]=unique(Prot_list);
USeq=List_seq(index_p,:);
Uniprot_id=Data(:,5);
Uni_prot=Uniprot_id(index_p);
clear index

%==========================================================================

%% Drug Similarities
%%1-Finger print (CS Sim)
[Feature1,NOT_find]=chem_feature_api(Drugs');
Feature_drug_finger=double(Feature1(:,33:end-7)); %%remove the 4 byte lenght prefix at first and 7 bit at end
[~,D_finger_simil]=Jaccard_sim_cal(Drugs,Feature_drug_finger); %      1- Jaccrad Similarity
S_D{1}=D_finger_simil;
%%2-Feature Based (FB Sim) 
[Feature_drug_ch]=All_drug_chemical_feature(Drugs);
[~,D_Feat_simil]=chemical_sim_cal(Drugs,Feature_drug_ch);
S_D{2}=D_Feat_simil;
%%3-TF-IDF (TF-IDF Sim) 
% [~,Feature_drug_TF5,D_TF_simil5]=TFidf_similarity(USmile',Drugs,5);
f_name1=strcat(data_name,'_D_TF_simil');
f_name3=strcat(f_name1,num2str(5));
%save(f_name3,'D_TF_simil5');
load(f_name3);
S_D{3}=D_TF_simil5;
%%4-LINGO (LINGO Sim)
% % % [D_Lingo_Sim5]=LingO_similarity(Drugs,Feature_drug_TF5);
f_name2=strcat(data_name,'_D_Lingo_Sim');
f_name4=strcat(f_name2,num2str(5));
%save(f_name3,'D_Lingo_Sim5');
load(f_name4);
S_D{4}=D_Lingo_Sim5;
for i=1:4
    [x,y]=find(isnan(S_D{i}));
    temp=S_D{i};
    for j=1:numel(x)
        temp(x(j),y(j))=0;
    end
    S_D{i}=temp;
end

%==========================================================================

%% Protein Similiarity
%%1-SW (SW Sim)
[~,P_SW_simil]=Swaling_sim_cal(Prots,USeq,1);
S_P{1}=P_SW_simil;
%%2-Feature Based (FB Sim)
Feature_Prot_ch=ALL_Protein_features(Prots,USeq,'my_set');
for i=1:size(Feature_Prot_ch,2)
    m1=max(Feature_Prot_ch(:,i));
    m2=min(Feature_Prot_ch(:,i));
    Feature_Prot_ch(:,i)=(Feature_Prot_ch(:,i)-m2)./(m1-m2);
end
[~,P_Feat_simil]=chemical_sim_cal(Prots,Feature_Prot_ch);
S_P{2}=P_Feat_simil;
%%3-TF-IDF (TF-IDF SIM)
% % [~,Feature_Prot_TF5,P_TF_simil5]=TFidf_similarity(USeq',Prots,5); 
f_name1=strcat(data_name,'_P_TF_simil');
f_name3=strcat(f_name1,num2str(5));
load(f_name3);
S_P{3}=P_TF_simil5;
%%4-LINGO (LINGO Sim)
% % [P_Lingo_Sim5]=LingO_similarity(Prots,Feature_Prot_TF5);
f_name2=strcat(data_name,'_P_Lingo_Sim');
f_name4=strcat(f_name2,num2str(5));
load(f_name4);
S_P{4}=P_Lingo_Sim5;
for i=1:4
    [x,y]=find(isnan(S_P{i}));
    temp=S_P{i};
    for j=1:numel(x)
        temp(x(j),y(j))=0;
    end
   S_P{i}=temp;     
end

%--------------------------------------------------------------------------

%% Indexes
for i=1:size(Data,1)
    ind_d(i)=find(ismember(Drugs,Drug_list(i))~=0);
    ind_p(i)=find(ismember(Prots,Prot_list(i))~=0);
end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [list,Similarity]=chemical_sim_cal(list,F)
[list,ind]=unique(list);
F=F(ind,:);
d=pdist2(F,F);
Similarity=1./(1+d);

end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [list,Similarity]=Jaccard_sim_cal(list,F)
[list,ind]=unique(list);
F=F(ind,:);
Similarity=1-pdist2(F,F,'jaccard');
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

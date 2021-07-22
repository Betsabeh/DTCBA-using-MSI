function [Interactions1]=cal_interaction_sim(Y,KD,ind_d,ind_p)
% This function calculate the Interaction Profile (IP)Similarity
%--------------------------------------------------------------------------
%[Interactions1]=cal_interaction_sim(Y,KD,ind_d,ind_p)
% Ipnuts:
% Y:        Interaction label according to the binarization threshold(labels +1,-1)
% KD:       Binding affinity values
% ind_d:    Drug indexes related to the drug-drug similarities
% ind_p:    Protein indexes related to the drug-drug similarities
% Output:
% Interactions1 : D*T Inteaction matrix 
%--------------------------------------------------------------------------
global Prots Drugs S_D S_P
num_P=numel(Prots);
num_D=numel(Drugs);
Interactions=zeros(num_D,num_P);
Interactions1=ones(num_D,num_P)*-5;%100000;
for i=1:numel(ind_d)
    Interactions(ind_d(i),ind_p(i))=Y(i);
    Interactions1(ind_d(i),ind_p(i))=KD(i);
end
%--------------------------------------------------------------------------
d=pdist2(Interactions',Interactions','Jaccard');
P_Inter_simil= exp(-d);
[x,y]=find(isnan(P_Inter_simil));
for i=1:numel(x)
    P_Inter_simil(x(i),y(i))=0;
end
S_P{2}=P_Inter_simil;

d=pdist2(Interactions,Interactions,'Jaccard');
D_Inter_simil = exp(-d);
[x,y]=find(isnan(D_Inter_simil));
for i=1:numel(x)
    D_Inter_simil(x(i),y(i))=0;
end
S_D{2}=D_Inter_simil;


end
%---------------------------------------------------------------------------------

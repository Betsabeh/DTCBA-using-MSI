function [Feature]=All_drug_chemical_feature(P_cid)
% This fuction finds the Features of Drugs
%---------------------------------------------------------------------------
% Download the Drug features from https://pubchem.ncbi.nlm.nih.gov/  according to CID 
% These features are in xlsx format 
%order of data in xlsx file :
%cid	cmpdname	cmpdsynonym	mw	mf	polararea	complexity	xlogp	heavycnt	hbonddonor	hbondacc	rotbonds	inchikey	iupacname	meshheadings	annothits	annothitcnt	aids	cidcdate
% selected index:
%mw=4 polararea=6 complexity=7 heavycnt=9 hbonddonor=10	hbondacc=11
%rotbonds=12 cidcdate=19
%-------------------------------------------------------------------------
[num,txt,raw]=xlsread('ALL_drug_chemical_properties.xlsx'); 
Pubchem_CID=num(:,1);
atoms={'C';'H';'N';'O';'F';'S';'Cl';'Br';'I'};%'Q';'Na';'P';'Br';'Si';'I';'B'};
n=numel(atoms);
k=1;
for i=1:numel(P_cid)
    id=find(Pubchem_CID==str2num(P_cid{i}));
    if isempty(id)
        NOT_found{k}=P_cid{i};
        k=k+1;
        continue;
    end
    Feature(i,1)=num(id,4);
    Feature(i,2)=num(id,6);
    Feature(i,3)=num(id,7);
    Feature(i,4)=num(id,9);
    Feature(i,5)=num(id,10);
    Feature(i,6)=num(id,11);
    Feature(i,7)=num(id,12);
% % % % %    t{i}=cell2mat(txt(id+1,5));
    Feature(i,8:7+n)=count_atoms(cell2mat(txt(id+1,5)));%my_chemical_find(cell2mat(txt(id+1,5)));%
 %   Feature(i,8)=num(id,19);
end
% NOT_found
end
%--------------------------------------------------------------------------
function F=count_atoms(chformula)
[num_atoms, atoms, species] = atomic(chformula);
atoms_array={'C';'H';'N';'O';'F';'S';'Cl';'Br';'I'};%'Q';'Na';;'P';'Br';'Si';'I';'B'};
F=zeros(1,size(atoms_array,1));
for i=1:numel(atoms)
    ind=find(ismember(atoms_array,atoms{i})~=0);
    if isempty(ind)==0
       F(ind)=num_atoms(i);
    else
        1
    end
        
end
end

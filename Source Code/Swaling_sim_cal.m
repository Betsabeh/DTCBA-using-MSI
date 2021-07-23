function [list,Similarity]=Swaling_sim_cal(list,F,type)
global data_name
switch type
    case 1
[list,ind_p]=unique(list);
file_name=strcat(data_name,'_Bind_Prot_Seq_SW_Sim.txt');
tt=textread(file_name,'%s','delimiter','\n','bufsize',30000);
Pr=cell2mat(tt(1,:));
Pr=strsplit(Pr,'		');
tt=tt(2:end,:);
P_simil=str2num(cell2mat(tt));
clear tt

for i=1:numel(list)
    ind_p(i)=find(ismember(Pr,list{i})~=0);
end
Similarity=P_simil(ind_p,ind_p);

%########################################################################
    case 2

file_name=strcat(data_name,'_Bind_Prot_Seq_SW_Sim.txt');
[list,ind]=unique(list);
F=F(ind);
for i=1:numel(F)
    t1=swalign(F{i},F{i});
    for j=i:numel(F)
        j
        t2=swalign(F{j},F{j});
        Similarity(i,j)=swalign(F{i},F{j})/(sqrt(t1)*sqrt(t2));
        Similarity(j,i)=Similarity(i,j);
    end
end
 write_similarity(file_name,list,Similarity)
end
end
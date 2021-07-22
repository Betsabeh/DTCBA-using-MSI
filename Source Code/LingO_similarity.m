function [Similarity]=LingO_similarity(P_Cid,TF)
%% calculate LINGO similarity based on TF-IDF
[UPcid,index,~]=unique(P_Cid);
d=size(TF,2);
Similarity=0;
for i=1:d
t1=pdist2(TF(:,i),TF(:,i),'cityblock');
t2=pdist2(TF(:,i),-1*TF(:,i),'cityblock');
temp=(t1./t2);
[x,y]=find(isnan(temp)==1);
for j=1:numel(x)
    temp(x(j),y(j))=0;
end
Similarity=Similarity+(1-temp);
end
Similarity=Similarity./d;
end



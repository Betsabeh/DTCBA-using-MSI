function [Uid,TF_ALL,Sim]=TFidf_similarity(Seq,id,q)
%This function calculates the TF-IDF similarity
%% q=qchar lingo
[Uid,index,~]=unique(id);
U_Seq=Seq(:,index);
[Voc,TF,Idf]=Vocab(U_Seq,q);
Sim=cal_sim(TF,Idf);
for i=1:numel(id)
    ind=find(ismember(Uid,id{i})~=0);
   TF_ALL(i,:)=TF(ind,:);
end

end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [Voc,TF,Idf]=Vocab(Seq,q)
k=1;
Voc=[];
for i=1:size(Seq,2)
    Doc=Seq{i};
    Doc=strrep(Doc,' ','');
    lenght=numel(Doc);
    TF(i,1)=0;
    for j=1:lenght-(q-1)
        temp=Doc(j:j+(q-1));
        ind=find(ismember(Voc,temp)~=0);
        if isempty(ind)==1
            Voc{k}=temp;
            TF(i,k)=1;
            k=k+1;
        else
            TF(i,ind)=TF(i,ind)+1;
        end
        %Doc_Token{i,j}=temp;
    end
end
N=size(TF,1);
% % % % % % % % % % % % % % % % % % % % hist(sum(TF))
df=TF~=0;
df=sum(df);

% % % % % ind=find(df<N/2);
% % % % % df=df(ind);
% % % % % TF=TF(:,ind);
% % % % % 
% % % % % % % % % % ind=find(df>3);
% % % % % % % % % % df=df(ind);
% % % % % % % % % % TF=TF(:,ind);

Idf=log2(N./df);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function Siml=cal_sim(TF,Idf)
num_Doc=size(TF,1);
for i=1:num_Doc
    m=max(TF(i,:));
    TF_Idf(i,:)=(TF(i,:)./m).*(Idf);
end
for i=1:size(TF_Idf,1)
    TF_Idf(i,:)=TF_Idf(i,:)./norm(TF_Idf(i,:));
end
Siml=1-pdist2(TF_Idf,TF_Idf,'cosine');
end
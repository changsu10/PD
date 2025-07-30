function Closet_Distance_ZScore_split(Cancer_Type,alpha,start,endn); %% MZ changed

load Data_mat/Gene_Distance  %%%%% distance between genes
load Data_mat/Net_PPI

%%  ---- MZ changed: moved data loading out of loop ----
load Data_mat/Map_List 
folder_path="./";
file_type=".txt";
final_file=strcat(folder_path,Cancer_Type,file_type);%input genes file name
gene_1 = load(final_file); %%%AK

load(['Data_mat/Gene_Drug_10uM']);
load Data_mat/Map_List
%% end of change

Net=Net(:,1:2);
Net(find(Net(:,1)-Net(:,2)==0),:)=[];
[LG,L]=largest_component(Net);
Net=LG{find(L==max(L))};
Genes=unique(Net(:));
Degree=Genes;
for i=1:length(Degree)
    Degree(i,2)=length(find(Net(:)==Degree(i)));
end
tic;
if nargin < 2
    alpha=0.5;
end

save_f = [num2str(start),'_network_proximity.txt'];%% MZ changed: save intermediate results
f = fopen(save_f, 'w');%% MZ changed

for II=start:endn %% MZ changed: only loop a subset to run paralell
    [z,p]=S_AB_Cal_Block(II,Net,Degree,Genes,Distance,gene_1,Map_List,Drug_List,Gene_Drug); %% MZ changed
    drug=Map_List{II,1};
    fprintf(f, '%s\t%f\t%f\n',drug, z, p);%% MZ changed: save intermediate results while running
    toc
    disp(II)
end
fclose(f);%% MZ changed

function [Z_Score,P_Value]=S_AB_Cal_Block(II,Net,Degree,Genes,Distance,gene_1,Map_List,Drug_List,Gene_Drug);%% MZ changed

drug=Map_List{II,1};
[i1,i2]=ismember(drug,Drug_List);
ii=i2;
gene_2=Gene_Drug(find(Gene_Drug(:,2)==ii),1);  %%%% drug gene
gene_1=intersect(gene_1,Genes);
gene_2=intersect(gene_2,Genes);
if isempty(gene_1) | isempty(gene_2)
    Z_Score=nan;
    P_Value=nan;%% MZ changed
else
    for i=1:1000
        g_1=drug_network_random_module_calculation(Degree,gene_1,i);
        g_2=drug_network_random_module_calculation(Degree,gene_2,i);
        CD(i,1)=CD_Cal(Genes,Distance,g_1,g_2);
    end
    CD(i+1,1)=CD_Cal(Genes,Distance,gene_1,gene_2);  %%%%% the first column
    Z_Score=(CD(end)-mean(CD(1:end-1)))/std(CD(1:end-1));
    [mu,sigma]=normfit(CD(1:end-1));%% MZ changed: get ref distribution mean and std
    P_Value=normcdf(CD(end),mu,sigma); %% MZ changed: calculated one-sided p-value
end

function CD=CD_Cal(Gene,Distance,gene_1,gene_2);
[a,m1]=ismember(gene_1,Gene);
[a,m2]=ismember(gene_2,Gene);  %%% drug gene
d=Distance(m1,m2);
d=double(d);
CD=mean([min(d)]);

function Module_R=drug_network_random_module_calculation(Degree,Gene,II);
ctime=datestr(now,30);
tseed=str2num(ctime((end-5):end));
rand('seed',tseed*II);    %%%%% random seed reset
Degree=sortrows(Degree,-2);
Degree_Dis=unique(Degree(:,2));
for i=1:length(Degree_Dis)
    Degree_Dis(i,2)=length(find(Degree(:,2)==Degree_Dis(i,1)));
end
Degree_Dis=sortrows(Degree_Dis,-1);
DB=Degree_Dis(1,1)+1;  %%% degree bin
kk=Degree_Dis(1,2);
for i=2:length(Degree_Dis)
    kk=kk+Degree_Dis(i,2);
    if kk>200
        DB=[DB;Degree_Dis(i,1)];
        kk=0;
    end
end
for i=1:length(DB)-1
    x=find(Degree(:,2)<DB(i) & Degree(:,2)>=DB(i+1));
    NDB{i,1}=Degree(x,1);   %%%% node in each degree bin
end
Module_R=[];
[a,b]=ismember(Gene,Degree(:,1));
MBD=Degree(b(a),2);
MBD=[Gene,MBD];
for i=1:length(Gene)
    k=MBD(i,2);
    m=find(DB<=k,1,'first')-1;
    NB=NDB{m};
    Module_R=[Module_R;NB(ceil(rand*length(NB)))];
end

%%%%%%% calculate the largest component of the PPI network
function [LG,L]=largest_component(G);
k=0;
if isempty(G)
    LG={[]};
    L=0;
end
while ~isempty(G)
    k=k+1;
    nn=G(ceil(rand*size(G,1)),1);
    m1=ismember(G(:,1),nn);
    m2=ismember(G(:,2),nn);
    g=G(m1|m2,:);
    G(m1|m2,:)=[];
    ng=unique(g(:));
    while ~isempty(ng)
        m1=ismember(G(:,1),ng);
        m2=ismember(G(:,2),ng);
        gg=G(m1|m2,:);
        ng=unique(gg(:));
        g=[g;gg];
        G(m1|m2,:)=[];
    end
    LG{k,1}=g;
    L(k,1)=length(unique(g(:)));
    clear g
end
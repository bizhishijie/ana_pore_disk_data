% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
load('p.mat');p=p';p_length=length(p);
for nn=1:30
    %% 加载数据
    disp(nn)
    load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
    load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
    load(['..\basic\pack\' fileList(nn).name '\edge.mat'])
    load(['..\basic\pack\' fileList(nn).name '\voro.mat'],'v','pore_rc')
    load(['..\basic\pack\' fileList(nn).name '\pore_1.mat'],'path_cell','is_del','neighbor_cell','ks');
    %%
    p=zeros(length(pore_rc),length(pore_rc));% 权重矩阵
    for ii=1:length(pore_rc)
        if ~is_del(ii)
            for jj=1:length(neighbor_cell{ii})
                pt=v(path_cell{ii}{jj},:);% path temp
                p(ii,neighbor_cell{ii}(jj))=sum(vecnorm(diff(pt)'));
                p(neighbor_cell{ii}(jj),ii)=sum(vecnorm(diff(pt)'));
            end
        end
    end
    p(p==0)=inf;
    G = graph(p);
    P = shortestpath(G,20,300);
    pore_tmp=20;
    path_tmp={};
   for ii=2:length(P)
        path_tmp{end+1}=path_cell{pore_tmp}(neighbor_cell{pore_tmp}==P(ii));
        pore_tmp=P(ii);
   end
   path_tmp(1)=[];% 删掉第一个
   for ii=1:length(path_tmp)
       x=v(cell2mat(path_tmp{ii}),1);
       y=v(cell2mat(path_tmp{ii}),2);
       z=v(cell2mat(path_tmp{ii}),3);
       plot3(x,y,z)
       hold on
   end
end
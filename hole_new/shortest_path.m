% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
load('p_intersect.mat');p=p';p_length=length(p);% 这里加载用于计算相交的采样点
for nn=1:30
    %% 加载数据
    disp(nn)
    load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
    load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
    load(['..\basic\pack\' fileList(nn).name '\edge.mat'])
    load(['..\basic\pack\' fileList(nn).name '\voro.mat'],'v','pore_rc')
    load(['..\basic\pack\' fileList(nn).name '\pore_1.mat'],'path_cell','is_del','neighbor_cell','ks','pore_c_id');
    %%
    power=inf(length(pore_rc),length(pore_rc));% 权重矩阵
    volume=nan(1,length(pore_rc));
    parfor ii=1:length(pore_rc)
        if ~is_del(ii)
            volume(ii)=volume_pore(ks{ii},pore_c_id{ii},v,Rc,Ori,p);
        end
    end
    pore_is_in=cellfun(@(c)all(is_in(c)),pore_c_id);
    volume_new=volume(pore_is_in'&~is_del);
    %     hist(volume_new)
    save(['..\basic\pack\' fileList(nn).name '\volume.mat'],'volume');
    for ii=1:length(pore_rc)
        if ~is_del(ii)
            for jj=1:length(neighbor_cell{ii})
                pt=v(path_cell{ii}{jj},:);% path temp
                power(ii,neighbor_cell{ii}(jj))=sum(vecnorm(diff(pt)'));
                power(neighbor_cell{ii}(jj),ii)=sum(vecnorm(diff(pt)'));
            end
        end
    end
    G = graph(power);
    %% 设定起点与终点并计算路线
    node_start=500;
    node_end=300;
    [P,d] = shortestpath(G,node_start,node_end);
    %     pore_tmp=node_start;
    %     path_tmp={};
    %     for ii=2:length(P)
    %         path_tmp{end+1}=path_cell{pore_tmp}(neighbor_cell{pore_tmp}==P(ii));
    %         pore_tmp=P(ii);
    %     end
    %     path_tmp(1)=[];% 删掉第一个

    %     for ii=1:length(path_tmp)
    %         x=v(cell2mat(path_tmp{ii}),1);
    %         y=v(cell2mat(path_tmp{ii}),2);
    %         z=v(cell2mat(path_tmp{ii}),3);
    %         plot3(x,y,z)
    %         hold on
    %     end
end
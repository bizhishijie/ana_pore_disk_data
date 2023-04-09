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
    power=zeros(length(pore_rc),length(pore_rc));% 权重矩阵
    volu=nan(1,length(pore_rc));
    parfor ii=1:length(pore_rc)
        warning('off')
        if ~is_del(ii)
            volu(ii)=volume_pore(ks{ii},pore_c_id{ii},v,Rc,Ori,p);
        end
    end
    save(['..\basic\pack\' fileList(nn).name '\volume.mat'],'volu');
end
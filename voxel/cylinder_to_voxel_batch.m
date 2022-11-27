% 将圆盘化为体素
clear  
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
voxel_num=1000;

parfor ii=1:length(fileList)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Ori').Ori;
    is_in=load(['..\basic\pack\' fileList(ii).name '\rc_is_in.mat']).is_in;
    disp(fileList(ii).name)

    Rc_min=min(Rc,[],2);
    Rc_max=max(Rc,[],2);% 体素划分的区域
    
    Rc=Rc(:,is_in);
    Ori=Ori(:,is_in);

    voxel=cylinder_to_voxel(Rc,Ori,[Rc_min Rc_max],cylinder_size,voxel_num);
    parsave(['..\basic\pack\' fileList(ii).name '\voxel.mat'],voxel)
end
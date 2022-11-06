% 将圆盘化为体素
clear  
d=30.75*2;
h=4.81*2;
r=d/2;
cylinder_size=[r,h];
fileList=dir('..\basic\pack\pack*');
voxel_num=500;

for ii=1:lengt(fileList)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Ori').Ori;
    disp(fileList(ii).name)

    Rc_min=min(Rc,[],2)+r;
    Rc_max=max(Rc,[],2)-r;

    voxel=cylinder_to_voxel(Rc,Ori,[Rc_min Rc_max],cylinder_size,voxel_num);
    save(['..\basic\pack\' fileList(ii).name '\voxel.mat'],'voxel')
end
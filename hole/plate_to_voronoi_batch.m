% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
load('tri.mat')
parfor ii=1:length(fileList)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Ori').Ori;
    disp(fileList(ii).name)
    [vertice,regions]=plate_to_voronoi(Rc,Ori,p);
    parsave(['..\basic\pack\' fileList(ii).name '\voronoi.mat'], vertice,regions);
end
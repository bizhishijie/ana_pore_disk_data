% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
d=30.75*2;
h=4.81*2*1.1;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%

parfor ii=1:length(fileList)
    disp(fileList(ii).name)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat']).Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat']).Ori;
    figure(ii)
    show_cylinder(Rc,Ori)

    Rc_max_sphere=max(Rc,[],2)-2*r;
    Rc_min_sphere=min(Rc,[],2)+2*r;

    Rc_2=Rc(:,all((Rc<Rc_max_sphere)&(Rc>Rc_min_sphere)));
    dt=delaunayTriangulation(Rc_2');
    polygon=convexHull(dt);% 质心组成的多边形

    k=convexHull(dt);

    trisurf(k, dt.Points(:,1),dt.Points(:,2),dt.Points(:,3), 'FaceColor', 'cyan','FaceAlpha',0.5)
    saveas(gcf,['..\basic\pack\' fileList(ii).name '\basic.fig'])
    saveas(gcf,['..\basic\pack\' fileList(ii).name '\basic.jpg'])
end
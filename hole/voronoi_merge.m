% 合并一个圆盘上的采样点的格子
clear
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];
fileList=dir('..\basic\pack\pack*');
load('tri.mat')
p_length=size(p,2);
for ii=2:length(fileList)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Ori').Ori;
    %     edge=load(['..\basic\pack\' fileList(ii).name '\edge.mat'],'edge').edge;
    %     rc_is_in=load(['..\basic\pack\' fileList(ii).name '\rc_is_in.mat'],'is_in').is_in;
    disp(fileList(ii).name)

    [vertice,regions]=cylinder_to_voronoi(Rc,Ori,p);
    point_cell=cell(1,length(Rc));
    volume_list=zeros(1,length(Rc));
    parfor jj=1:length(Rc)
        disp(jj);
        % 对每一个圆盘遍历，对应regions中的前若干个点
        id_tmp=regions((jj-1)*p_length+1:jj*p_length);

        id_tmp=cell2mat(cellfun(@transpose,id_tmp,'UniformOutput',false));% 获得全部的id列表
        id_tmp=sort(id_tmp);
        id_uni=unique(id_tmp);
        id_statistics=tabulate(id_tmp);
        id_statistics=id_statistics(:,1:2);
        id_statistics(id_statistics(:,2)==0,:)=[];
        id_uni(id_statistics(:,2)>=4)=[];
        % 在 id  _tmp 里面，出现4次及以上的是内部的点，应该删除
        cylinder_point=vertice(id_uni,:);
        %             [V_inter,cylinder_point]=intersects_polyhedron_(cylinder_point',edge);
        % 找到对应的点的编号
        if ~isempty(cylinder_point)&&~any(any(isinf(cylinder_point)))
            [k,volume] = convhull(cylinder_point,"Simplify",true);
            dt_tmp=delaunayTriangulation(cylinder_point);
            point_cell{jj}=dt_tmp.Points;
            volume_list(jj)=volume;
        else
            dt_tmp=[];
            volume=0;
            point_cell{jj}=[];
            volume_list(jj)=0;
        end
        %     draw_cylinder(Rc(:,jj),Ori(:,jj),'r',r,h,0.5)
        %     trisurf(k,dt_tmp.Points(:,1),dt_tmp.Points(:,2),dt_tmp.Points(:,3),'FaceAlpha',0.5)

    end
    parsave(['..\basic\pack\' fileList(ii).name '\voronoi_point.mat'],point_cell);
    parsave(['..\basic\pack\' fileList(ii).name '\voronoi_volume.mat'],volume_list);
end
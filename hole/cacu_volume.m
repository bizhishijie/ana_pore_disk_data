clear
fileList=dir('..\basic\pack\pack*');
load('tri.mat')
p_length=size(p,2);

d=30.75*2;
h=4.81*2;
r=d/2;
cylinder_size=[r,h];
parfor ii=1:length(fileList)
    disp(fileList(ii).name)
    cylinder_point_cell=load(['..\basic\pack\' fileList(ii).name '\voronoi_cut_cylinder_point.mat']).cylinder_point_cell;
    %     edge_cell=load(['..\basic\pack\' fileList(ii).name '\voronoi_cut_edge.mat']).edge_cell;
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat']).Rc;
    Rc_min=min(Rc,[],2)+r;
    Rc_max=max(Rc,[],2)-r;
    box=[Rc_min Rc_max];
    volume_list=zeros(1,length(cylinder_point_cell));
    for jj=1:length(cylinder_point_cell)
        if ~isempty(cylinder_point_cell{jj})
            if is_in_box(cylinder_point_cell{jj},box)
                tri = delaunayTriangulation(cylinder_point_cell{jj});
                [ch ,volume] = convexHull(tri);
                volume_list(jj)=volume;
            end
        end
    end
    parsave(['..\basic\pack\' fileList(ii).name '\voronoi_volume.mat'],volume_list)
end
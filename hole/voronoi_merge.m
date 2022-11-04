clear
fileList=dir('..\basic\pack\pack*');
load('tri.mat')
p_length=size(p,2);
parfor ii=1:length(fileList)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Ori').Ori;
    disp(fileList(ii).name)
    %     [vertice,regions]=plate_to_voronoi(Rc,Ori,p);
    %     parsave(['..\basic\pack\' fileList(ii).name '\voronoi.mat'], vertice,regions);
    regions=load(['..\basic\pack\' fileList(ii).name '\voronoi.mat'],'regions').regions;
    vertice=load(['..\basic\pack\' fileList(ii).name '\voronoi.mat'],'vertice').vertice;
    edge_cell=cell(1,size(Rc,2));
    cylinder_point_cell=cell(1,size(Rc,2));
    for jj=1:size(Rc,2)
        % 对每一个圆盘遍历，对应regions中的前若干个点
        id_tmp=regions((jj-1)*p_length+1:jj*p_length);
        id_tmp=cell2mat(cellfun(@transpose,id_tmp,'UniformOutput',false));
        id_tmp=sort(id_tmp);
        id_uni=unique(id_tmp);
        id_statistics=zeros(length(id_uni),2);
        id_statistics(:,1)=id_uni;
        for kk=1:size(id_uni)
            id_statistics(kk,2)=sum(id_tmp==id_uni(kk));
        end
        id_uni(id_statistics(:,2)>=4)=[];
        % 在 id  _tmp 里面，出现4次及以上的是内部的点，应该删除
        cylinder_point=vertice(id_uni,:);
        % 找到对应的点的编号
        if ~any(isinf(cylinder_point))
            edge = convhull(cylinder_point);
        else
            continue
        end
        %     defaultFaceColor  = [0.6875 0.8750 0.8984];
        %     trisurf(edge, cylinder_point(:,1),cylinder_point(:,2),cylinder_point(:,3) , ...
        %         'FaceColor', defaultFaceColor, 'FaceAlpha',0.8)
        %     hold on
        %     disp(jj)
        edge_cell{jj}=edge;
        cylinder_point_cell{jj}=cylinder_point;
    end
    parsave(['..\basic\pack\' fileList(ii).name '\voronoi_cut_cylinder_point.mat'], cylinder_point_cell);
    parsave(['..\basic\pack\' fileList(ii).name '\voronoi_cut_edge.mat'], edge_cell);
end
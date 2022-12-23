% 读取圆柱的质心坐标，并输出voronoi的参数
% 并合并一个圆盘上的采样点的格子
clear
d=30.75*2/0.8;
h=4.81*2/0.8;% 避免穿模
r=d/2;% 改这里没有用，请更改生成p.mat的文件
cylinder_size=[r,h];
fileList=dir('..\basic\pack\pack*'); 
load('p.mat');p=p';
p_length=size(p,2);
for ii=1:length(fileList)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Ori').Ori;
    %     edge=load(['..\basic\pack\' fileList(ii).name '\edge.mat'],'edge').edge;
    %     rc_is_in=load(['..\basic\pack\' fileList(ii).name '\rc_is_in.mat'],'is_in').is_in;
    disp(fileList(ii).name)

    [vertice,regions]=cylinder_to_voronoi(Rc,Ori,p);
    point_cell=cell(1,length(Rc));
    parfor jj=1:length(Rc)
        disp(jj);
        % 对每一个圆盘遍历，对应regions中的前若干个点
        id_tmp=regions((jj-1)*p_length+1:jj*p_length);
        
%         hold on
%         for nn=1:p_length
%             dt=delaunayTriangulation(vertice(  cell2mat(regions((jj-1)*p_length+nn)) , :  ));
%             tetramesh(dt,'FaceColor','cyan');
%         end

        id_tmp=cell2mat(cellfun(@transpose,id_tmp,'UniformOutput',false));% 获得全部的id列表
        id_tmp=sort(id_tmp);
        id_uni=unique(id_tmp);
        id_statistics=tabulate(id_tmp);
        id_statistics=id_statistics(:,1:2);
        id_statistics(id_statistics(:,2)==0,:)=[];
        id_uni(id_statistics(:,2)>=4)=[];
        % 在 id  _tmp 里面，出现4次及以上的是内部的点，应该删除
        cylinder_point=vertice(id_uni,:);% 取出属于这个圆盘的边界点
        % 找到对应的点的编号
        if ~isempty(cylinder_point)&&~any(any(isinf(cylinder_point)))
            cylinder_point=unique(cylinder_point,'rows');% 由于相邻的边界必有相同点，这里作unique
            dt=delaunayTriangulation(cylinder_point);% 不可以只取出最外侧的点
            point_cell{jj}=dt.Points;
        else
            point_cell{jj}=[];
        end
        %     draw_cylinder(Rc(:,jj),Ori(:,jj),'r',r,h,0.5)
        %     trisurf(k,dt_tmp.Points(:,1),dt_tmp.Points(:,2),dt_tmp.Points(:,3),'FaceAlpha',0.5)

    end
    save(['..\basic\pack\' fileList(ii).name '\voronoi_point.mat'],'point_cell','vertice','regions');
end
% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
nn=1;

load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
load(['..\basic\pack\' fileList(nn).name '\p_all.mat'])
load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
load(['..\basic\pack\' fileList(nn).name '\point_pore.mat'])% 每一个孔的位置
load(['..\basic\pack\' fileList(nn).name '\voronoi_point.mat'])
load('p.mat');p=p';
disp(fileList(nn).name);
Rc=Rc(:,is_in);Ori=Ori(:,is_in);% 由于 near_pore_cylinder_id 的编号是筛选后的编号，这里需要提前筛选
p_all=p_all(is_in);
point_cell=point_cell(is_in);
% near_pore_cylinder 是距离每一个 pore 最近的圆柱的编号
pore_neighbor_cell=cell(1,length(near_pore_cylinder_id));
[near_pore_cylinder_id,ia,ic]=unique(near_pore_cylinder_id','rows');% 删去接触相同圆盘的无效孔
near_pore_cylinder_id=near_pore_cylinder_id';% 每一个圆柱的邻居
point_pore=point_pore(:,ia);
for ii=1:length(point_pore)
    pore(ii).rc=point_pore(:,ii);
    pore(ii).near_cylinder=near_pore_cylinder_id(:,ii);
    pore(ii).neighbor=find(sum(ismember(near_pore_cylinder_id,near_pore_cylinder_id(:,ii)))==3);% 该点的邻居的id
    edge_unique=unique(cell2mat(cellfun(@transpose,point_cell(pore(ii).near_cylinder),'UniformOutput',false))','rows');
    pore(ii).voronoi_edge=edge_unique';
end
for ii=1:length(point_pore)
    cnt=0;
    real_neighbor=[];
    path=cell(0);
    for jj=pore(ii).neighbor
        cnt=cnt+1;
        common_cylinder=intersect(pore(ii).near_cylinder,pore(jj).near_cylinder);
        p1=cell2mat(cellfun(@transpose,point_cell(common_cylinder(1)),'UniformOutput',false));
        p2=cell2mat(cellfun(@transpose,point_cell(common_cylinder(2)),'UniformOutput',false));
        p3=cell2mat(cellfun(@transpose,point_cell(common_cylinder(3)),'UniformOutput',false));
        p_tmp=intersect(p1',p2','rows')';
        p_common=intersect(p_tmp',p3','rows')';
        if ~isempty(p_common) && ismember(pore(jj).rc',p_common','rows') && ismember(pore(ii).rc',p_common','rows')
            real_neighbor(end+1)=jj;
            path{end+1}=p_common;
        end
    end
    pore(ii).real_neighbor=real_neighbor;
    pore(ii).path=path;
end
show_pore(pore,2,Rc,Ori)
hold on;clf
% plot3(p1(1,:),p1(2,:),p1(3,:),'r.')
% plot3(p2(1,:),p2(2,:),p2(3,:),'g.')
% plot3(p3(1,:),p3(2,:),p3(3,:),'b.')
plot3(p_common(1,:),p_common(2,:),p_common(3,:),'m.')


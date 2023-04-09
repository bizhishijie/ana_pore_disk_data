% 根据划分的voronoi计算孔的位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];
fileList=dir('..\basic\pack\pack*');
for nn=1:length(fileList)
    Rc=load(['..\basic\pack\' fileList(nn).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(nn).name '\basic.mat'],'Ori').Ori;
    point_cell=load(['..\basic\pack\' fileList(nn).name '\voronoi_point.mat'],'point_cell').point_cell;
    edge=load(['..\basic\pack\' fileList(nn).name '\edge.mat'],'edge').edge;
    is_in=load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'],'is_in').is_in;
    disp(fileList(nn).name)
    Rc=Rc(:,is_in);Ori=Ori(:,is_in);
    point_cell=point_cell(is_in);

    point_list=cell2mat(cellfun(@transpose,point_cell,'UniformOutput',false));
    idx_face=convhulln(edge');
    rc2_c=mean(edge,2);
    is_in_list=false(1,length(point_list));
    parfor ii=1:length(point_list)
        face_side=zeros(1,size(idx_face,1));
        for jj=1:size(idx_face,1)
            n_vec=cross(edge(:,idx_face(jj,1))-edge(:,idx_face(jj,2)),edge(:,idx_face(jj,2))-edge(:,idx_face(jj,3)));
            n_vec=n_vec/norm(n_vec);
            face_side_tmp=dot(point_list(:,ii)-edge(:,idx_face(jj,1)),n_vec)*dot(rc2_c-edge(:,idx_face(jj,1)),n_vec);
            if face_side_tmp>0
                face_side(jj)=1;
            else
                break
            end
        end
        is_in_list(ii)=(sum(face_side==0)==0);
    end
    % 删去在外部的voronoi边界点

    is_in_list=find(is_in_list);
    point_list=point_list(:,is_in_list)';
    [point_list_unique,m,n]=unique(point_list,'rows');
    point_statistics=tabulate(n);
    point_statistics=sortrows(point_statistics,2);
    point_pore=point_list_unique(point_statistics(point_statistics(:,2)==4,1),:)';% 统计频次为4的点
    near_pore_cylinder_id=cell(1,length(point_pore));

    parfor ii=1:length(point_pore)
        near_pore_cylinder_id{ii}=[];
        for jj=1:length(Rc)
            if ismember(point_pore(:,ii)',point_cell{jj},"rows")
                near_pore_cylinder_id{ii}(end+1)=jj;
            end
        end
        near_pore_cylinder_id{ii}=near_pore_cylinder_id{ii}';
    end% 找到距离每个孔的圆盘
    near_pore_cylinder_id=cell2mat(near_pore_cylinder_id);
    save(['..\basic\pack\' fileList(nn).name '\point_pore.mat'],'point_pore','near_pore_cylinder_id');
    % 编号为通过is_in筛选后的编号
end
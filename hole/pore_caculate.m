% 根据划分的voronoi计算孔的位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];
fileList=dir('..\basic\pack\pack*');
load('p.mat')
p_length=size(p,2);
for nn=1:length(fileList)
    Rc=load(['..\basic\pack\' fileList(nn).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(nn).name '\basic.mat'],'Ori').Ori;
    point_cell=load(['..\basic\pack\' fileList(nn).name '\voronoi_point.mat'],'point_cell').point_cell;
    edge=load(['..\basic\pack\' fileList(nn).name '\edge.mat'],'edge').edge;
    is_in=load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'],'is_in').is_in;
    disp(fileList(nn).name)

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
        % 删去在外部的voronoi边界点
        is_in_list(ii)=(sum(face_side==0)==0);
    end
    
    is_in_list=find(is_in_list);
    point_list=point_list(:,is_in_list)';
    [point_list_unique,m,n]=unique(point_list,'rows');
    point_statistics=tabulate(n);
    point_statistics=sortrows(point_statistics,2);
    point_pore=point_list_unique(point_statistics(point_statistics(:,2)==4,1),:)';% 统计频次为4的点

    %     Rc=Rc';
    %     plot3(Rc(is_in,1),Rc(is_in,2),Rc(is_in,3),'bo');
    %     hold on
    %     plot3(point_pore(:,1),point_pore(:,2),point_pore(:,3),'ro');
    %     axis equal

    Rc=Rc(:,is_in);Ori=Ori(:,is_in);
    d=zeros(length(Rc),length(point_pore));
    for ii=1:length(Rc)
        rc=Rc(:,ii);ori=Ori(:,ii);
        Ac=[rc+ori*h/2 rc-ori*h/2];
        parfor jj=1:length(point_pore)
            [Distance,~,~] = distance_cylinder_sphere(Ac, r, point_pore(:,jj), 0);
            d(ii,jj)=Distance;% 计算每一个孔的点到各个圆盘的距离
        end
    end

    [~,near_pore_cylinder_id]=mink(d,4);% 计算距离最近的4个点
    save(['..\basic\pack\' fileList(nn).name '\point_pore.mat'],'point_pore','d','near_pore_cylinder_id');
    % 编号为通过is_in筛选后的编号
end
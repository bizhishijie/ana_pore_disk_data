clear
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
for nn=1:length(fileList)
    load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
    load(['..\basic\pack\' fileList(nn).name '\p_all.mat'])
    load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
    load(['..\basic\pack\' fileList(nn).name '\point_pore.mat'])% 每一个孔的位置
    load('p.mat');
    disp(fileList(nn).name);
    Rc=Rc(:,is_in);Ori=Ori(:,is_in);% 由于编号都是筛选后的编号，这里需要提前筛选
    p_all=p_all(is_in);
    % near_pore_cylinder 是距离每一个 pore 最近的圆柱的编号
    pore_neighbor_cell=cell(1,length(near_pore_cylinder_id));
    [near_pore_cylinder_id,ia,ic]=unique(near_pore_cylinder_id','rows');% 删去接触相同圆盘的无效孔
    near_pore_cylinder_id=near_pore_cylinder_id';% 每一个圆柱的邻居
    point_pore=point_pore(:,ia);
    pore_area_cell=cell(1,length(near_pore_cylinder_id));
    parfor ii =1:length(near_pore_cylinder_id)
        pore_neighbor=find(sum(ismember(near_pore_cylinder_id,near_pore_cylinder_id(:,ii)))==3);% 该点的邻居的id
        % 和该孔有相同的三个相邻圆盘的孔的id，定义为该孔的邻居
        pore_area=zeros(12,3*length(pore_neighbor));
        % 第一个维度的1~3、4~6、7~9分别表示中垂面的两个向量，10~12表示中垂面经过的点
        pore_neighbor_cell{ii}=pore_neighbor;
        for jj=1:length(pore_neighbor)% 对于每一个该孔的邻居计算
            A1=[point_pore(:,pore_neighbor(jj)) point_pore(:,ii)];
            % 从孔到其邻居的向量
            common_cylinder=intersect(near_pore_cylinder_id(:,pore_neighbor(jj)), ...
                near_pore_cylinder_id(:,ii));
            % 和他们共同的圆柱的编号
            %         for kk=1:length(common_cylinder)
            for kk=1:3
                rc=Rc(:,common_cylinder(kk));
                ori=Ori(:,common_cylinder(kk));
                A0=[rc-ori*h/2 rc+ori*h/2];% 共同的圆柱的从底到顶的向量
                [Distance,P0,P1] = distance_cylinders(A0, r, A1, 0);% 计算得到的P0是圆柱上的点，P1是向量上的点
                plane=null((A1(:,2)-A1(:,1))');% 两孔的中垂面的两个向量
                pore_area(1:6,(jj-1)*3+kk)=plane(:);
                pore_area(7:9,(jj-1)*3+kk)=(A1(:,2)-A1(:,1))'/norm(A1(:,2)-A1(:,1));
                pore_area(10:12,(jj-1)*3+kk)=P1;
                % P1是连线最近的点
            end
        end
        pore_area_cell{ii}=pore_area;
    end
    % sum(isempty(pore_neighbor)); % 没有空的
    % 计算各个嗓子所围成的形状
    cube=[0 1 1;0 1 -1; 0 -1 1; 0 -1 -1;
        2 1 1;2 1 -1; 2 -1 1; 2 -1 -1];
    edge_cell=cell(1,length(pore_area_cell));

    normal_vec=[0;0;1];
    length_p=size(p,1);
    for ii = 1:length(pore_area_cell)
        pore_area=pore_area_cell{ii};
        edge=(cube)'*10*r;
        for jj=1:size(pore_area,2)
            base=reshape(pore_area(1:9,jj),3,3);
            cube_tmp=base*edge;
            k = convhull(edge','Simplify',true);
            trisurf(k,edge(1,:),edge(2,:),edge(3,:),'FaceColor','cyan','FaceAlpha',0.5)

            [~,edge1]=intersects_polyhedron_(edge,cube_tmp+pore_area(10:12,jj)+pore_area(7:9,jj)*r);
            % 新生成的和原来的求交集
            if ~isempty(edge1)
                edge=edge1;
            end
        end
        pore_neighbor=pore_neighbor_cell{ii};
        for jj=1:size(near_pore_cylinder_id,1)% 对于每一个该孔的相邻的圆柱计算
            p_all_cylinder=p_all{near_pore_cylinder_id(jj)};% 加载该圆柱的采样点
            for kk=1:length(p_all)
                if inpolygon3d_new(p_all(:,kk),edge)
                    edge=[edge p_all(:,kk)];
                end
            end
        end
        edge_cell{ii}=edge;
    end
    save(['..\basic\pack\' fileList(nn).name '\pore_area.mat'],'pore_area_cell','near_pore_cylinder_id','edge_cell')% 并保存剩余的孔的id
end
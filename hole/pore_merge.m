% 读取圆盘位置
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
    load(['..\basic\pack\' fileList(nn).name '\voronoi_point.mat'])
    load('p.mat');p=p';p_length=length(p);
    disp(fileList(nn).name);
    Rc=Rc(:,is_in);Ori=Ori(:,is_in);% 由于 near_pore_cylinder_id 的编号是筛选后的编号，这里需要提前筛选
    p_all=p_all(is_in);
    point_cell=point_cell(is_in);
    % near_pore_cylinder 是距离每一个 pore 最近的圆柱的编号
    %%
    % 计算孔的属性
    pore_neighbor_cell=cell(1,length(near_pore_cylinder_id));
    [near_pore_cylinder_id,ia,ic]=unique(near_pore_cylinder_id','rows');% 删去接触相同圆盘的无效孔
    near_pore_cylinder_id=near_pore_cylinder_id';% 每一个圆柱的邻居
    point_pore=point_pore(:,ia);
    for ii=1:length(point_pore)
        pore(ii).rc=point_pore(:,ii);
        dis=distance_point_cylinder(point_pore(:,ii),Rc,Ori,cylinder_size);
        pore(ii).rk=min(dis(1,:));% 距离最近的圆柱的距离
        pore(ii).near_cylinder=near_pore_cylinder_id(:,ii);
        pore(ii).neighbor=find(sum(ismember(near_pore_cylinder_id,near_pore_cylinder_id(:,ii)))==3);% 该点的邻居的id
        pore(ii).is_del=false;
    end
    %%
    % 计算孔到孔的路径
    for ii=1:length(point_pore)
        cnt=0;
        path={};
        for jj=pore(ii).neighbor
            cnt=cnt+1;
            common_cylinder=intersect(pore(ii).near_cylinder,pore(jj).near_cylinder);
            p1=cell2mat(cellfun(@transpose,point_cell(common_cylinder(1)),'UniformOutput',false));
            p2=cell2mat(cellfun(@transpose,point_cell(common_cylinder(2)),'UniformOutput',false));
            p3=cell2mat(cellfun(@transpose,point_cell(common_cylinder(3)),'UniformOutput',false));
            % 三个圆盘的voronoi边界的点
            p_common=intersect(intersect(p1',p2','rows'),p3','rows')';
            if ~isempty(p_common) && ismember(pore(jj).rc',p_common','rows') && ismember(pore(ii).rc',p_common','rows')
                path{end+1}=p_common;
            end
        end
        pore(ii).path=path;
    end
    %%
    for ii=1:length(pore)
        if pore(ii).is_del
            continue
        end
        pore(ii).throat=zeros(length(pore(ii).path),6);% 分别保存嗓子的位置和方向
        for jj=1:length(pore(ii).path)
            path_tmp=pore(ii).path{jj};
            dis=inf(1,length(path_tmp));
            for kk=1:size(path_tmp,2)
                if ~is_point_in_cylinder(path_tmp(:,kk),Rc,Ori,cylinder_size)
                    dis_tmp=distance_point_cylinder(path_tmp(:,kk),Rc,Ori,cylinder_size);
                    dis(kk)=min(dis_tmp(1,:));
                end
            end
            [~,idx]=min(dis);% 找嗓子的位置
            pore(ii).throat(jj,1:3)=path_tmp(:,idx)';
            if idx==1
                pore(ii).throat(jj,4:6)=path_tmp(:,idx+1)'-path_tmp(:,idx)';
            else
                pore(ii).throat(jj,4:6)=path_tmp(:,idx)'-path_tmp(:,idx-1)';
            end
        end
    end
    %%
    map=1:length(is_in);
    map=map(is_in);
    for ii=1:length(point_pore)
        disp(ii)
        P={};% 保存每个孔的采样点的凸包
        K={};
        neighbor=pore(ii).neighbor;
        for nei=1:length(neighbor)% 对每一个颗粒的每一个邻居遍历
            jj=neighbor(nn);
            common_cylinder=intersect(pore(ii).near_cylinder,pore(jj).near_cylinder);
            for c1=1:length(common_cylinder)% 遍历圆柱的组合
                cylinder_id_1=map(common_cylinder(c1));% 圆柱的id
                rc_id_1=common_cylinder(c1);% 对应的质心的id

                id_1=regions((cylinder_id_1-1)*p_length+1:cylinder_id_1*p_length);
                id_1_list=cell2mat(cellfun(@transpose,id_1,'UniformOutput',false));
                p_tmp_1=cellfun(@(P)vertice(P,:),id_1,'UniformOutput',false);% 采样点的voronoi边界的点
                try
                    con_tmp_1=cellfun(@(P)convhulln(vertice(P,:)),id_1,'UniformOutput',false);% 采样点的voronoi边界的点组成的凸包
                catch
                    continue
                end
                p_1=cylinder_to_point(Rc(:,rc_id_1),Ori(:,rc_id_1),p);% 圆柱的采样点
                for c2=c1+1:length(common_cylinder)
                    cylinder_id_2=map(common_cylinder(c2));% 圆柱的id
                    rc_id_2=common_cylinder(c2);% 对应的质心的id
                    id_2=regions((cylinder_id_2-1)*p_length+1:cylinder_id_2*p_length);
                    id_2_list=cell2mat(cellfun(@transpose,id_2,'UniformOutput',false));
                    p_tmp_2=cellfun(@(P)vertice(P,:),id_2,'UniformOutput',false);% 采样点的voronoi边界的点
                    try

                        con_tmp_2=cellfun(@(P)convhulln(vertice(P,:)),id_2,'UniformOutput',false);
                    catch
                        continue
                    end
                    p_2=cylinder_to_point(Rc(:,rc_id_2),Ori(:,rc_id_2),p);% 圆柱的采样点
                    id_common=intersect(id_1_list,id_2_list);% 共同的点的id，只需要找包含这个编号的采样点的voronoi
                    if isempty(id_common)
                        continue% 如果没有相同的点则继续
                    end

                    f=@(id_cell)any(ismember(id_common,id_cell));% 哪个voronoi需要展示出来
                    is_show_1=cell2mat(cellfun(f,id_1,'UniformOutput',false));
                    is_show_2=cell2mat(cellfun(f,id_2,'UniformOutput',false));
                    p_1_show=p_tmp_1(is_show_1);
                    p_2_show=p_tmp_2(is_show_2);

                    for cc=1:length(p_1_show)
                        p_1_show_tmp=unique(p_1_show{cc},'rows');
                        flag=true(length(p_1_show_tmp),1);
                        for tt=1:size(pore(ii).throat,1)
                            throat=pore(ii).throat(tt,:);
                            flat=@(x,y,z)(throat(4)*(x-throat(1))+throat(5)*(y-throat(2))+throat(6)*(z-throat(3)));
                            flag_standard=sign(flat(pore(ii).rc(1),pore(ii).rc(2),pore(ii).rc(3)));% 孔的质心带入throat的方程
                            is_same_side=sign(flat(p_1_show_tmp(:,1),p_1_show_tmp(:,2),p_1_show_tmp(:,3)))==flag_standard;
                            flag=flag&is_same_side;
                        end
                        p_1_show_tmp(~flag,:)=[];
                        try
                            common_1_tmp=convhulln(p_1_show_tmp);
                            K{end+1}=common_1_tmp;P{end+1}=p_1_show_tmp;
                        catch
                            continue
                        end
                    end
                    for cc=1:length(p_2_show)
                        p_2_show_tmp=unique(p_2_show{cc},'rows');
                        flag=true(length(p_2_show_tmp),1);
                        for tt=1:size(pore(ii).throat,1)
                            throat=pore(ii).throat(tt,:);
                            flat=@(x,y,z)(throat(4)*(x-throat(1))+throat(5)*(y-throat(2))+throat(6)*(z-throat(3)));
                            flag_standard=sign(flat(pore(ii).rc(1),pore(ii).rc(2),pore(ii).rc(3)));% 孔的质心带入throat的方程
                            is_same_side=sign(flat(p_2_show_tmp(:,1),p_2_show_tmp(:,2),p_2_show_tmp(:,3)))==flag_standard;
                            flag=flag&is_same_side;
                        end
                        p_2_show_tmp(~flag,:)=[];
                        try
                            common_2_tmp=convhulln(p_2_show_tmp);
                            K{end+1}=common_2_tmp;P{end+1}=p_2_show_tmp;
                        catch
                            continue
                        end
                    end
                end
            end
        end
        pore(ii).connective=K;
        pore(ii).area=P;
    end
    %%
    % 合并
    flag=true;
    while flag
        disp(1)
        flag=false;
        for ii=1:length(pore)
            if pore(ii).is_del
                continue
            end
            cnt=0;% 一共已经计算了几个pore
            del_num=0;% 该pore中已经删去的元素数
            for jj=pore(ii).neighbor
                cnt=cnt+1;
                if pore(jj).is_del
                    continue
                end
                if norm(pore(ii).rc-pore(jj).rc)<=pore(ii).rk && pore(jj).rk<=pore(ii).rk% 如果满足合并条件
                    flag=true;
                    pore(jj).is_del=true;% jj被合并
                    neighbor0=pore(ii).neighbor;
                    pore(ii).connective=[pore(ii).connective pore(jj).connective];
                    pore(ii).area=[pore(ii).area pore(jj).area];
                    neighbor_new=setdiff(pore(jj).neighbor,pore(ii).neighbor);% 新的邻居
                    neighbor_new([pore(neighbor_new).is_del])=[];% 已经被删除的邻居就不算了
                    neighbor_new(neighbor_new==ii)=[];% 自己不是自己的邻居
                    pore(ii).neighbor(cnt-del_num)=[];% 删除jj
                    pore(ii).neighbor=[pore(ii).neighbor neighbor_new];% 更新新的邻居，放到最后
                    %                     pore(ii).voronoi_edge=unique([pore(ii).voronoi_edge pore(jj).voronoi_edge]','rows')';% 更新新的边界点
                    pore(ii).near_cylinder=unique([pore(ii).near_cylinder;pore(jj).near_cylinder]);% 更新新的相邻的圆柱

                    path0=pore(ii).path{cnt-del_num};% 原本的 path
                    pore(ii).path(cnt-del_num)=[];% 删除原本的path
                    for kk=1:length(pore(jj).neighbor)
                        if ~pore(pore(jj).neighbor(kk)).is_del && pore(jj).neighbor(kk)~=ii && ~ismember(pore(jj).neighbor(kk),neighbor0)
                            % 如果路径到达的pore没被删除并且路径的终不是ii并且邻居是新的
                            pore(ii).path{end+1}=[path0 pore(jj).path{kk}];
                        end
                    end% 更新新的路径
                    del_num=del_num+1;
                end
            end
        end
    end
    %%
    % 计算喉咙
    for ii=1:length(pore)
        if pore(ii).is_del
            continue
        end
        pore(ii).throat=zeros(length(pore(ii).path),6);% 分别保存嗓子的位置和方向
        for jj=1:length(pore(ii).path)
            path_tmp=pore(ii).path{jj};
            dis=inf(1,length(path_tmp));
            for kk=1:size(path_tmp,2)
                if ~is_point_in_cylinder(path_tmp(:,kk),Rc,Ori,cylinder_size)
                    dis_tmp=distance_point_cylinder(path_tmp(:,kk),Rc,Ori,cylinder_size);
                    dis(kk)=min(dis_tmp(1,:));
                end
            end
            [~,idx]=min(dis);% 找嗓子的位置
            pore(ii).throat(jj,1:3)=path_tmp(:,idx)';
            if idx==1
                pore(ii).throat(jj,4:6)=path_tmp(:,idx+1)'-path_tmp(:,idx)';
            else
                pore(ii).throat(jj,4:6)=path_tmp(:,idx)'-path_tmp(:,idx-1)';
            end
        end
    end
    save(['..\basic\pack\' fileList(nn).name '\pore.mat'],'pore')
end
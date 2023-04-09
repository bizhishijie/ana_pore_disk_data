load('D:\毕个业\basic\pack\pack104\pore.mat')
load('D:\毕个业\basic\pack\pack104\rc_is_in.mat')
load('D:\毕个业\basic\pack\pack104\basic.mat')
load('D:\毕个业\basic\pack\pack104\voronoi_point.mat')
Rc=Rc(:,is_in);Ori=Ori(:,is_in);
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];%
load('p.mat');p=p';
p_length=length(p);
%%
% 计算每个邻居所对应的圆柱的具有相同点的采样点的voronoi边界的并集减去圆柱
% is_in 前后的圆柱编号对应
map=1:length(is_in);
map=map(is_in);
voronoi_cell=cell(1,length(pore));
for ii=1:length(pore)
    if pore(ii).is_del
        continue
    end
    disp(ii)
    neighbor=pore(ii).neighbor;
    P={};% 保存每个孔的采样点的凸包
    K={};
    for nn=1:length(neighbor)% 对每一个颗粒的每一个邻居遍历
        jj=neighbor(nn);
        common_cylinder=intersect(pore(ii).near_cylinder,pore(jj).near_cylinder);
        for c1=1:length(common_cylinder)% 遍历圆柱的组合
            cylinder_id_1=map(common_cylinder(c1));% 圆柱的id
            rc_id_1=common_cylinder(c1);% 对应的质心的id
            id_1=regions((cylinder_id_1-1)*p_length+1:cylinder_id_1*p_length);
            id_1_list=cell2mat(cellfun(@transpose,id_1,'UniformOutput',false));
            p_tmp_1=cellfun(@(P)vertice(P,:),id_1,'UniformOutput',false);% 采样点的voronoi边界的点
            con_tmp_1=cellfun(@(P)convhulln(vertice(P,:)),id_1,'UniformOutput',false);% 采样点的voronoi边界的点组成的凸包

            p_1=cylinder_to_point(Rc(:,rc_id_1),Ori(:,rc_id_1),p);% 圆柱的采样点
            for c2=c1+1:length(common_cylinder)
                cylinder_id_2=map(common_cylinder(c2));% 圆柱的id
                rc_id_2=common_cylinder(c2);% 对应的质心的id
                id_2=regions((cylinder_id_2-1)*p_length+1:cylinder_id_2*p_length);
                id_2_list=cell2mat(cellfun(@transpose,id_2,'UniformOutput',false));
                p_tmp_2=cellfun(@(P)vertice(P,:),id_2,'UniformOutput',false);% 采样点的voronoi边界的点
                con_tmp_2=cellfun(@(P)convhulln(vertice(P,:)),id_2,'UniformOutput',false);

                p_2=cylinder_to_point(Rc(:,rc_id_2),Ori(:,rc_id_2),p);% 圆柱的采样点
                id_common=intersect(id_1_list,id_2_list);% 共同的点的id，只需要找包含这个编号的采样点的voronoi
                if isempty(id_common)
                    continue% 如果没有相同的点则继续
                end
                %                 point_list_common=vertice(id_common,:);%共同的点可以画出来看
                f=@(id_cell)any(ismember(id_common,id_cell));% 哪个voronoi需要展示出来
                is_show_1=cell2mat(cellfun(f,id_1,'UniformOutput',false));
                is_show_2=cell2mat(cellfun(f,id_2,'UniformOutput',false));
                p_1_show=p_tmp_1(is_show_1);
                p_2_show=p_tmp_2(is_show_2);
                % 分别找到对应的voronoi采样点的id

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
    pore(ii).point_show=P;
end
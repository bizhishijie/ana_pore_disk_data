clear;clc
load('D:\毕个业\basic\pack\pack12\voro.mat')
load('D:\毕个业\basic\pack\pack12\pore_1.mat')
load('D:\毕个业\basic\pack\pack12\volume.mat')
load('D:\毕个业\basic\pack\pack12\basic.mat')
load('D:\毕个业\basic\pack\pack12\edge.mat')
load('D:\毕个业\basic\pack\pack12\rc_is_in.mat')
d1=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d1/2;
cylinder_size=[r,h];%
%%
volu(volu<0)=0;
pore_is_in=cellfun(@(c)all(is_in(c)),pore_c_id);
volume_new=volu(pore_is_in'&~is_del);
% hist(volume_new)
range=[max(edge,[],2),min(edge,[],2)];
cnt=0;
% path_id=cell2mat(cellfun(@(c)cell2mat(c'),path_cell,'UniformOutput',false)');
power=inf(length(pore_rc),length(pore_rc));% 权重矩阵
for ii=1:length(pore_rc)
    if ~is_del(ii)
        for jj=1:length(neighbor_cell{ii})
            pt=v(path_cell{ii}{jj},:);% path temp
            dl=sum(vecnorm(diff(pt)'));
            power(ii,neighbor_cell{ii}(jj))=dl;
            power(neighbor_cell{ii}(jj),ii)=dl;
        end
    end
end
G = graph(power);
%%
simu_num=1e3;
rate_list=zeros(1,simu_num);
while cnt<simu_num% 共模拟多少个孔中的点
    %% 对每次采样遍历
    l=0;% 路径的总长度
    p1=rand(3,1).*(range(:,1)-range(:,2))+range(:,2);
    p2=rand(3,1).*(range(:,1)-range(:,2))+range(:,2);
    if inpolygon3d_new(p1,edge) && ~is_point_in_cylinder(p1,Rc,Ori,cylinder_size) && inpolygon3d_new(p2,edge) && ~is_point_in_cylinder(p2,Rc,Ori,cylinder_size)
        cnt=cnt+1;
        % 随机到的点在edge内，且不在圆柱内
        % 找到距离每个点最近的path点
        d1=distance_point_cylinder(p1,Rc,Ori,cylinder_size);
        [~,idx_c1]=mink(d1(1,:),5);% 取最近的若干圆柱看他们的采样点
        idx_p=cellfun(@(n)any(ismember(idx_c1,n)),pore_c_id);
        p1_near=v(cell2mat([path_cell{idx_p}]'),:);
        d2=distance_point_cylinder(p2,Rc,Ori,cylinder_size);
        [~,idx_c2]=mink(d2(1,:),5);% 取最近的若干圆柱看他们的采样点
        idx_p=cellfun(@(n)any(ismember(idx_c2,n)),pore_c_id);
        p2_near=v(cell2mat([path_cell{idx_p}]'),:);
        [d_1_1,p1_n_id]=min(vecnorm(p1_near'-p1));
        [d_1_2,p2_n_id]=min(vecnorm(p2_near'-p2));% 找到距离随机的点的最近的路径上的点
        %% 画图
        %         show_cylinder(Rc(:,idx_c1),Ori(:,idx_c1))
        %         show_cylinder(Rc(:,idx_c2),Ori(:,idx_c2))
        %         plot3(p1(1),p1(2),p1(3),'o');plot3(p2(1),p2(2),p2(3),'o');
        %         plot3(p1_near(:,1),p1_near(:,2),p1_near(:,3),'o');plot3(p2_near(:,1),p2_near(:,2),p2_near(:,3),'o');
        %%
        l=l+d_1_1+d_1_2;% 这部分比较近，忽略掉优化
        ps_1=find(sum(p1_near(p1_n_id,:)==v,2)==3);
        ps_2=find(sum(p2_near(p2_n_id,:)==v,2)==3);
        is_1=false(1,length(path_cell));
        is_2=false(1,length(path_cell));
        parfor ii=1:length(path_cell)
            if any(cellfun(@(c)(ismember(ps_1,c)),path_cell{ii}))
                is_1(ii)=true;
            end
            if any(cellfun(@(c)(ismember(ps_2,c)),path_cell{ii}))
                is_2(ii)=true;
            end
        end
        is_1=find(is_1);is_2=find(is_2);
        d=inf(length(is_1),length(is_2));
        for ii=1:length(is_1)
            for jj=1:length(is_2)
                [~,d(ii,jj)] = shortestpath(G,is_1(ii),is_2(jj));
            end
        end
        if all(isinf(d),'all')
            cnt=cnt-1;
            continue
        end
        l=l+min(d,[],'all');
        [ii,jj]=find(d==min(d,[],'all'));
        path_tmp_1=path_cell{is_1(ii)}{cellfun(@(c)(ismember(ps_1,c)),path_cell{is_1(ii)})};
        l=l+sum(vecnorm(diff(v(path_tmp_1(1:find(path_tmp_1==ps_1)),:))'));
        path_tmp_2=path_cell{is_2(jj)}{cellfun(@(c)(ismember(ps_2,c)),path_cell{is_2(jj)})};
        l=l+sum(vecnorm(diff(v(path_tmp_2(1:find(path_tmp_2==ps_2)),:))'));
        rate_list(cnt)=l/norm(p1-p2);
    end
end
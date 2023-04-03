% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
% load('p.mat');p=p';p_length=length(p);
for nn=1:30
    %% 加载数据
    disp(nn)
    load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
    load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
    load(['..\basic\pack\' fileList(nn).name '\edge.mat'])
    load(['..\basic\pack\' fileList(nn).name '\voro.mat'],'v','pore_rc')
    load(['..\basic\pack\' fileList(nn).name '\volume.mat'])
    load(['..\basic\pack\' fileList(nn).name '\pore_1.mat'],'ks','pore_c_id','is_del')
    load(['..\basic\pack\' fileList(nn).name '\idx_wrong_disk.mat'])
    %%
    volu(volu<0)=0;
    pore_is_in=cellfun(@(c)all(is_in(c)),pore_c_id);
    volume_new=volu(pore_is_in'&~is_del);
    range=[max(edge,[],2),min(edge,[],2)];
    %%
    simu_num=1e6;% 共模拟多少个孔中的点
    rate_list=zeros(1,simu_num);
    p_list=zeros(6,simu_num);
    near_idx=cell(1,simu_num);
    distance_list=zeros(1,simu_num);
    parfor pp=1:simu_num
        %% 对每次采样遍历
        warning('off')
        p1=rand(3,1).*(range(:,1)-range(:,2))+range(:,2);
        flag=~(inpolygon3d_new(p1,edge) && ~is_point_in_cylinder(p1,Rc,Ori,cylinder_size));
        while flag
            flag=~(inpolygon3d_new(p1,edge) && ~is_point_in_cylinder(p1,Rc,Ori,cylinder_size));
            p1=rand(3,1).*(range(:,1)-range(:,2))+range(:,2);
        end
        %         disp(pp)

        d1=distance_point_cylinder(p1,Rc,Ori,cylinder_size);
        [~,idx_c1]=mink(d1(1,:),4);% 取最近的若干圆柱看他们的采样点
        if ismember(idx_c1,idx_wrong_disk)
            continue
        end
        idx_p=find(cellfun(@(n)all(ismember(idx_c1,n)),pore_c_id)&~is_del');% 采样点附近的孔的编号
        is_in_pore=false(1,length(idx_p));
        for ii=1:length(idx_p)
            P=v(unique(ks{idx_p(ii)}),:);
            is_in_pore(ii)=inpolygon3d_new(p1,P');
        end
        if any(idx_p(is_in_pore))
            near_idx{pp}=idx_p(is_in_pore);
            d=distance_point_cylinder(p1,Rc,Ori,cylinder_size);
            distance_list(pp)=min(d(1,:));
        end
        %         show_pore(pore_rc(pore_id,:),ks{pore_id},v,Rc(:,pore_c_id{pore_id}),Ori(:,pore_c_id{pore_id}),pore_rc(neighbor_cell{pore_id},:),path_cell(pore_id),throat_cell{pore_id});
        %         plot3(p1(1),p1(2),p1(3))
        %         sum(is_in_pore)
    end
    %%
    rc_pore=cell(1,length(pore_rc));
    for ii=1:simu_num
        if distance_list(ii)~=0
            for jj=1:length(near_idx{ii})
                rc_pore{near_idx{ii}(jj)}=[rc_pore{near_idx{ii}(jj)} distance_list(ii)];
            end
        end
    end
    %%
    save(['..\basic\pack\' fileList(nn).name '\rc_pore.mat'],'rc_pore');
end
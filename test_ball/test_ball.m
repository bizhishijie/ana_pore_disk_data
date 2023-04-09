clear;clc
r=23.7/2;
load('p_sphere.mat');p=(p*r*0.5)';n=length(p);
load('basic_418.mat');
[~,min_idx]=mink(vecnorm(Rc-mean(Rc,2)),20);
Rc=Rc(:,min_idx);
Td_map=[2,3,4;1,3,4;1,2,4;1,2,3];
Li_map=[1,2;1,3;1,4;2,3;2,4;3,4];

ball_sample_point=cell2mat(cellfun(@(c)(p+c),num2cell(Rc,1),'UniformOutput',false))';
Td=delaunayn(ball_sample_point);
Td_fix=floor((Td-1)/n)+1;

Td_outside=~cellfun(@(c)length(unique(c))==1,num2cell(Td_fix,2));
Td=Td(Td_outside,:);
Td_fix=Td_fix(Td_outside,:);% 删掉在颗粒内部的四面体
f=@(t,line_tmp,Rc_tmp)(-min(vecnorm((t*(line_tmp(1:3)-line_tmp(4:6))+line_tmp(4:6))-Rc)-r));

rs_pore=zeros(length(Td),3);
dist_pore=zeros(length(Td),1);
Nd_tmp_id=cell(1,length(Td));
Td_inside=false(1,length(Td));
%%
parfor ii=1:length(Td)
    warning('off')
    Rc_eff_id=Td_fix(ii,:);
    tetra_tmp=ball_sample_point(Td(ii,:),:);% 三维的对应的是四面体
    %     plot(tri_tmp(:,1),tri_tmp(:,2),'r')
    % 四面体
    DT = delaunayTriangulation(tetra_tmp);
    rs_pore_list1=circumcenter(DT);% 四面体的外心
    if ~isempty(rs_pore_list1)% 四点共面找不到外心
        dist_pore_list1=min(vecnorm(rs_pore_list1'-Rc))-r;
        flag_tetra=is_point_in_tetra(rs_pore_list1,tetra_tmp(1,:),tetra_tmp(2,:),tetra_tmp(3,:),tetra_tmp(4,:));
        rs_pore_list1=rs_pore_list1(flag_tetra,:);
        dist_pore_list1=dist_pore_list1(flag_tetra);
    else
        dist_pore_list1=[];
    end

    rs_pore_list2=zeros(4,3);
    dist_pore_list2=zeros(4,1);
    is_in_tri=false(4,1);
    is_in_tetra=false(4,1);
    is_in_diff_disk=false(4,1);
    % 如果不在四面体的内部则需要计算四个面的外心的距离
    for jj=1:size(Td_map,1)
        tri_tmp_cir=tetra_tmp(Td_map(jj,:),:);
        C = circumcenter_3d(tri_tmp_cir(1,:),tri_tmp_cir(2,:),tri_tmp_cir(3,:));

        rs_pore_list2(jj,:)=C;
        dist_pore_list2(jj)=min(vecnorm(C'-Rc))-r;
        is_in_tri(jj)=is_point_in_tri(C,tri_tmp_cir(1,:),tri_tmp_cir(2,:),tri_tmp_cir(3,:));
        is_in_tetra(jj)=is_point_in_tetra(C,tetra_tmp(1,:),tetra_tmp(2,:),tetra_tmp(3,:),tetra_tmp(4,:));
        is_in_diff_disk(jj)=length(unique(Td_fix(ii,Td_map(jj,:))))>1;
    end
    dist_pore_list2=dist_pore_list2(is_in_tri&is_in_diff_disk&is_in_tetra);
    rs_pore_list2=rs_pore_list2(is_in_tri&is_in_diff_disk&is_in_tetra,:);

    rs_pore_list3=zeros(6,3);
    dist_pore_list3=zeros(6,1);
    for jj=1:size(Li_map,1)
        line_tmp=reshape(tetra_tmp(Li_map(jj,:),:)',[],1);
        [min_point_position,dist_pore_list3(jj)]=fminbnd(@(t)f(t,line_tmp,Rc),0,1);
        rs_pore_list3(jj,:)=min_point_position*(line_tmp(1:3)-line_tmp(4:6))+line_tmp(4:6);
        is_in_diff_disk(jj)=length(unique(Td_fix(ii,Li_map(jj,:))))>1;
    end
    dist_pore_list3=-dist_pore_list3(is_in_diff_disk);
    [dist_pore(ii),max_idx]=max([dist_pore_list1;dist_pore_list2;dist_pore_list3]);
    rs_pore_tmp=[rs_pore_list1;rs_pore_list2;rs_pore_list3];
    rs_pore(ii,:)=rs_pore_tmp(max_idx,:);
    Td_inside(ii)=dist_pore(ii)<0;

    % figure(1);clf;hold on;axis equal
    % tetramesh(DT)
    % plot3(rs_pore(ii,1),rs_pore(ii,2),rs_pore(ii,3),'*')
    % show_ball(Rc(:,unique(Rc_eff_id)))
    % disp(dist_pore(ii))
end
%%
% 计算准pore的邻居
% Nd_tmp_id=cellfun(@(c)find(sum(ismember(Td,c),2)==3),num2cell(Td,2),'UniformOutput',false);\
Td=Td(~Td_inside,:);
Td_fix=Td_fix(~Td_inside,:);
parfor ii=1:length(Td)
    Nd_tmp_id{ii}=find(sum(ismember(Td,Td(ii,:)),2)==3)';
end
%% 第一次合并
IDX_merge=cell(length(Td),1);
dist_pore_tmp2=cellfun(@(c)dist_pore(c),Nd_tmp_id,'UniformOutput',false);
parfor ii=1:size(Td,1)
    n_idx0=Nd_tmp_id{ii}';
    dist_tmp=dist_pore_tmp2{ii};
    n_idx1=n_idx0(dist_tmp==max(dist_tmp));% 和最大的合并
    n_idx1=n_idx1(1);
    n_idx2=n_idx0(length(unique(Td_fix(ii,:)))>3&all(sort(Td_fix(ii,:))==sort(Td_fix(n_idx0,:),2),2));% 和完全相同的颗粒接触且是pore，直接合并
    IDX_merge{ii}=unique([n_idx1;n_idx2;ii]);
end% 每个三角形合并的三角形的id
IDX_merge1 = cell_combine_(IDX_merge);
% plot(cellfun('length',IDX_merge1))
%% 第二次合并
IDX_merge2=IDX_merge1;
flag=true;
while flag
    idx_tmp=1:length(IDX_merge2);
    idx_tmp=repelem(idx_tmp,cellfun('length',IDX_merge2));
    IDX_convert(cell2mat(IDX_merge2))=idx_tmp;% 每个四面体中心merge到的编号
    IDX_merge_pore=cell(1,length(IDX_merge2));
    for ii=1:length(IDX_merge2)% 第二次合并
        IDX_merge2_tmp=IDX_merge2{ii};
        merge_particle=unique(Td_fix(IDX_merge2_tmp,:));% 这个pore对应的颗粒
        pore_neighbor=setdiff(unique(cell2mat(Nd_tmp_id(IDX_merge2_tmp))),IDX_merge2_tmp);% pore的邻居
        merge_neighbor=IDX_convert(pore_neighbor);% merge的邻居
        merge_neighbor_particle=cellfun(@(c)unique(Td_fix(c,:)),IDX_merge2(merge_neighbor),'UniformOutput',false);
        is_neighbor=cellfun(@(c)all(ismember(merge_particle,c)),merge_neighbor_particle);% 不能带来新的颗粒邻居
        %     any(is_neighbor),length(merge_particle)
        pore_neighbor=pore_neighbor(is_neighbor);
        merge_neighbor=merge_neighbor(is_neighbor);

        [c,~,ic]=unique(merge_neighbor);
        
        neighbor_dist=dist_pore(pore_neighbor);
        merge2_dist=zeros(1,length(c));
        for jj=1:length(c)
            merge2_dist(jj)=min(neighbor_dist(jj==ic));
        end
        [~,max_idx]=max(merge2_dist);
        IDX_merge_pore{ii}=[c(max_idx),ii];
    end
    flag=any(cellfun('length',IDX_merge_pore)>1);
    IDX_merge_pore_1 = cell_combine_(IDX_merge_pore);
    IDX_merge2=cellfun(@(c)cell2mat(IDX_merge2(c)),IDX_merge_pore_1,'UniformOutput',false);
end
% plot(cellfun('length',IDX_merge2))
%%
dx=r/10;
rc_g_cell=cell(1,length(IDX_merge2));
idx_dcell_cell=cell(1,length(IDX_merge2));
% IDX_merge2 Td_fix Td cylinder_sample_point
% save(['../basic/pack/' fileList(nn).name '/pore.mat'],'IDX_merge2','Td_fix','Td','Cv','cylinder_sample_point')
for ii=1:length(IDX_merge2)
    %% 划分网格
    % disp(ii)
    idx_tetra=IDX_merge2{ii};
    idx_dcell=Td(idx_tetra,:);
    idx_dcell=unique(idx_dcell(:));
    Rc_tmp=Rc(:,unique(Td_fix(idx_tetra,:)));
    %         if any(~is_in(idx_dcell))
    %             continue
    %         end
    Rc_d=ball_sample_point(idx_dcell,:)';
    range=[min(Rc_d,[],2) max(Rc_d,[],2)];
    [xg,yg,zg]=meshgrid(range(1,1):dx:range(1,2),range(2,1):dx:range(2,2),range(3,1):dx:range(3,2));
    rc_g=[xg(:) yg(:) zg(:)]';
    if length(rc_g)>1000000% 超出内存限制，多半是没有合并好
        continue
    end
    clear Rc_d idx_dcell xg yg zg
    %% 判断点是否需要留下
    idx_out_particle=(all(dist(rc_g',Rc)>r,2))';
    idx_in_tetra_loop=false(length(idx_tetra),size(rc_g,2));

    parfor jj=1:length(idx_tetra)
        Rc_tetra=ball_sample_point(Td(idx_tetra(jj),:),:)';
        A=Rc_tetra(:,1)-Rc_tetra(:,4);
        B=Rc_tetra(:,2)-Rc_tetra(:,4);
        C=Rc_tetra(:,3)-Rc_tetra(:,4);
        rc_g2=rc_g-Rc_tetra(:,4);

        M3=[A B C]'*[A B C];
        N3=[A B C]'*rc_g2;
        cf3=M3\N3;
        idx_in_tetra_loop(jj,:)=sum(cf3,1)>=0&sum(cf3,1)<=1&all(cf3>0,1)&all(cf3<1,1);
    end
    clear A B C M3 N3 cf3
    idx_in_tetra=any(idx_in_tetra_loop,1);
    idx_keep=idx_out_particle&idx_in_tetra;
    rc_g=rc_g(:,idx_keep);
    idx_dcell=unique(Td_fix(idx_tetra,:));

    clf
    show_sphere(Rc(:,idx_dcell))
    plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.');axis equal;hold on
    
    rc_g_cell{ii}=rc_g;
    idx_dcell_cell{ii}=idx_dcell;
end
cell_not_empty=~cellfun('isempty',rc_g_cell);
rc_g_cell=rc_g_cell(cell_not_empty);
idx_dcell_cell=idx_dcell_cell(cell_not_empty);
% cell_in=cellfun(@(c)all(~isempty(c)&c<range(:,1)&c>range(:,2),"all"),rc_g_cell);
cell_eff=cellfun('length',idx_dcell_cell)>3;
% rc_g_cell=rc_g_cell(cell_in&cell_eff);
rc_g_cell=rc_g_cell(cell_eff);
% idx_dcell_cell=idx_dcell_cell(cell_in&cell_eff);
idx_dcell_cell=idx_dcell_cell(cell_eff);
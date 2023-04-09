clear;clc;
load('p.mat')
d=30.75*2/0.8;
h=4.81*2/0.8;D=h;
r=d/2;
dx=h/20;
cylinder_size=[r,h];
% 读取圆盘位置
fileList=dir('..\basic\pack\pack*');
% for nn=1:30
nn=26;
%%
load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
load(['..\basic\pack\' fileList(nn).name '\edge.mat'])
load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
load(['..\basic\pack\' fileList(nn).name '\idx_wrong_disk.mat'])

Rc(:,idx_wrong_disk)=[];Ori(:,idx_wrong_disk)=[];
% Rc=Rc(:,1:500);Ori=Ori(:,1:500);is_in=is_in(1:500);
% plot3(p(:,1),p(:,2),p(:,3),'.')
Rc0=Rc;
% plot3(Rc0(1,:),Rc0(2,:),Rc0(3,:),'.')
Rc=cylinder_to_point(Rc,Ori,p');
[Vv,Cv]=voronoin(Rc');
% Td=DelaunayTri(Rc').Triangulation;
Td=delaunayn(Rc');
Td_fix=fix((Td-1)./length(p))+1;% 检查pore归属于某几个颗粒，+1是因为整除
% pore_fix_id=find(cellfun(@(c)length(unique(c)),num2cell(Td_fix,2))==4);% 在多个颗粒中间的id
% find(cellfun(@(c)length(unique(c)),num2cell(Td_fix,2))==1)% 不存在颗粒内部的
rs_pore=zeros(3,size(Td,1));
dist_pore=zeros(1,size(Td,1));
Nd_idx=cell(1,size(Td,1));
parfor ii=1:size(Td,1)
    freq=tabulate(cell2mat(Cv(Td(ii,:)).'));
    idx=freq(freq(:,2)==4,1);
    idx=idx(idx~=1);% to remove inf
    rs_pore_tmp=Vv(idx,:)';
    rs_pore(:,ii)=rs_pore_tmp;
    Rc_eff_id=unique(Td_fix(ii,:));
    dist_pore(ii)=min(distance_point_cylinder(rs_pore_tmp,Rc0(:,Rc_eff_id),Ori(:,Rc_eff_id),cylinder_size));
    %     Nd_idx{ii}=find(sum(ismember(Td,Td(ii,:)),2)==3);
    % 邻居的定义需被修改，不能是一个面的两边的两个点
    Nd_tmp_id=find(sum(ismember(Td,Td(ii,:)),2)==3);% possible neighbors
    Td_tmp=Td_fix(Nd_tmp_id,:);% possible neighbors' particle
    flag=sum(Td_tmp-mode(Td_fix(ii,:),2)==0,2)==3;% 看每一行最多次数的
    Nd_idx{ii}=Nd_tmp_id(~flag);% 判断不能跨行
    % 不能有三个共同的点在一个圆盘上且四个不同
end
%%
IDX_merge=cell(1,size(Td,1));
edge_pore=[min(Rc,[],2) max(Rc,[],2)];
edge_rate=0.1;
edge_pore=[edge_pore(:,1)+(edge_pore(:,2)-edge_pore(:,1))*edge_rate edge_pore(:,2)-(edge_pore(:,2)-edge_pore(:,1))*edge_rate];
pore_inside=find(all(rs_pore>edge_pore(:,1)&rs_pore<edge_pore(:,2)));
% pore_outside=find(any(rs_pore<edge_pore(:,1)|rs_pore>edge_pore(:,2)));% intersect(pore_inside,pore_outside)
%%
parfor ii=1:size(Td,1)% 向大的方向合并
    if any(ii==pore_inside)
        n_idx=Nd_idx{ii}';n_idx=intersect(n_idx,pore_inside);
        n_idx=n_idx(dist_pore(n_idx)==max(dist_pore(n_idx)));
        if length(n_idx)>1
            n_idx=n_idx(1);
        end
        IDX_merge{ii}=intersect(n_idx,pore_inside);
    end
end
IDX_merge1 = cell_combine_(IDX_merge);
idx_tmp=1:length(IDX_merge1);
idx_tmp=repelem(idx_tmp,cellfun('length',IDX_merge1));
IDX_convert(cell2mat(IDX_merge1))=idx_tmp;% 每个四面体中心merge到的编号
IDX_merge_pore=cell(1,length(IDX_merge1));
%%
parfor ii=1:length(IDX_merge1)% connect for the second time
    pore_neighbor=intersect(unique(cell2mat(Nd_idx(IDX_merge1{ii})')),pore_inside);% pore的邻居
    pore_particle=unique(Td_fix(IDX_merge1{ii},:));
    merge_neighbor=unique(IDX_convert(pore_neighbor));% merge的邻居
    merge_neighbor_particle=cellfun(@(c)unique(Td_fix(IDX_merge1{c},:)),num2cell(merge_neighbor),'UniformOutput',false);
    merge_neighbor=merge_neighbor(cellfun(@(c)all(ismember(c,pore_particle)),merge_neighbor_particle));
    merge_neighbor_dist=cellfun(@(c)max(dist_pore(IDX_merge1{c})),num2cell(merge_neighbor));
    merge_neighbor=merge_neighbor(merge_neighbor_dist==max(merge_neighbor_dist));
    if length(merge_neighbor)>1
        merge_neighbor=merge_neighbor(1);
    end
    IDX_merge_pore{ii}=merge_neighbor;
end
IDX_merge_pore_1 = cell_combine_(IDX_merge_pore);
IDX_merge2=cellfun(@(c)cell2mat(IDX_merge1(c)),IDX_merge_pore_1,'UniformOutput',false);
% 把在外面的删掉
%%
IDX_merge2=cellfun(@(c)intersect(c,pore_inside),IDX_merge2,'UniformOutput',false);
% IDX_merge2=cellfun(@(c)setdiff(c,pore_outside),IDX_merge2,'UniformOutput',false);
IDX_merge2=IDX_merge2(cell2mat(cellfun(@(c)length(unique(Td_fix(c,:)))>=4,IDX_merge2,'UniformOutput',false)));
IDX_merge2=IDX_merge2(~cellfun('isempty',IDX_merge2));
%%
% clf;plot(cellfun('length',IDX_merge2))
%%
% show connection
figure(1);clf;
ii=114;
idx_show=IDX_merge2{ii};
show_cylinder(Rc0(:,unique(Td_fix(idx_show,:))),Ori(:,unique(Td_fix(idx_show,:))),'r')
for ii=1:length(idx_show)
    if ~isempty(IDX_merge{idx_show(ii)})
        p_tmp=[rs_pore(:,IDX_merge{idx_show(ii)}) rs_pore(:,idx_show(ii))];
        plot3(p_tmp(1,:),p_tmp(2,:),p_tmp(3,:),'b-o')
        hold on
    else
        p_tmp=[rs_pore(:,idx_show(ii))];
        plot3(p_tmp(1,:),p_tmp(2,:),p_tmp(3,:),'o','MarkerFaceColor','r')
    end
end
xlabel('x');ylabel('y');zlabel('z')
axis equal
%%
clear pore_inside pore_outside
% for ii=1:length(IDX_merge2)
[xs,ys,zs]=sphere(30);
% h=surf(xs,ys,zs,ones(size(xs)));
r_sphere_cell=cell(1,length(IDX_merge2));
rc_g_cell=cell(1,length(IDX_merge2));
range_cell=cell(1,length(IDX_merge2));
for ii=1:length(IDX_merge2)
    idx_tetra=IDX_merge2{ii};
    idx_dcell=Td(idx_tetra,:);
    idx_dcell=unique(idx_dcell(:));
    Rc0_tmp=Rc0(:,unique(Td_fix(idx_tetra,:)));
    Ori_tmp=Ori(:,unique(Td_fix(idx_tetra,:)));
    %         if any(~is_in(idx_dcell))
    %             continue
    %         end
    Rc_d=Rc(:,idx_dcell);
    range=[min(Rc_d,[],2) max(Rc_d,[],2)];
    [xg,yg,zg]=meshgrid(range(1,1):dx:range(1,2),range(2,1):dx:range(2,2),range(3,1):dx:range(3,2));
    rc_g=[xg(:) yg(:) zg(:)]';
    if length(rc_g)>10000000% 超出内存限制，多半是没有合并好
        continue
    end
    clear Rc_d idx_dcell xg yg zg
    %%
    d1=dot((reshape(rc_g,3,1,[])-Rc0_tmp),repmat(Ori_tmp,1,1,size(rc_g,2)));
    d2=vecnorm((reshape(rc_g,3,1,[])-Rc0_tmp)-d1.*repmat(Ori_tmp,1,1,size(rc_g,2)));
    d1=max(abs(d1)-h/2,0);d2=max(d2-r,0);
    idx_out_particle=permute(all(d1>0|d2>0,[1,2]),[2,3,1]);% 在同一个圆盘的两个面的外面
    clear d1 d2
    idx_in_tetra_loop=zeros(length(idx_tetra),size(rc_g,2));
    parfor jj=1:length(idx_tetra)
        Rc_tetra=Rc(:,Td(idx_tetra(jj),:));
        A=Rc_tetra(:,1)-Rc_tetra(:,4);
        B=Rc_tetra(:,2)-Rc_tetra(:,4);
        C=Rc_tetra(:,3)-Rc_tetra(:,4);
        rc_g2=rc_g-Rc_tetra(:,4); 

        M3=[A B C]'*[A B C];
        N3=[A B C]'*rc_g2;
        cf3=M3\N3;
        idx_in_tetra_loop(jj,:)=sum(cf3,1)>=0&sum(cf3,1)<=1&all(cf3>0,1);
    end
    clear A B C M3 N3 cf3
    idx_in_tetra=any(idx_in_tetra_loop,1);
    idx_keep=idx_out_particle&idx_in_tetra;
    rc_g=rc_g(:,idx_keep);
    %         plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.')
    clear idx_out_particle idx_in_tetra idx_keep
    %%
    A=cell(length(idx_tetra));
    for aa=1:length(idx_tetra)
        for bb=1:length(idx_tetra)
            A{aa,bb}=intersect(Td(idx_tetra(aa),:),Td(idx_tetra(bb),:));
        end
    end
    surf_common=unique(cell2mat(A(cellfun(@length,A)==3)),'rows');
    d_tetra_loop=inf(4*length(idx_tetra),size(rc_g,2));
    for jj=1:length(idx_tetra)
        Rc_tetra=Rc(:,Td(idx_tetra(jj),:));
        neighbor_surf=Td(sum(ismember(Td,Td(idx_tetra(jj),:)'),2)==3,:);
        surf_common2=surf_common(all(sum(reshape(Td(idx_tetra(jj),:),1,1,[])==surf_common,3),2),:);
        % 四面体和它的邻居的共同的面，每行代表一个
        map=[];
        map(Td(idx_tetra(jj),:))=1:4;% 将面的编号挪到1~4
        idx_delete=mod(10-sum(map(surf_common2),2)+2,4);% the rows need to be deleted
        idx_delete(idx_delete==0)=4;
        idx_permutation=[1 2 3 4;2 3 4 1;3 4 1 2;4 1 2 3];
        idx_permutation(idx_delete,:)=[];% 删掉和它的邻居重复的面的组合
        for kk=1:size(idx_permutation,1)
            A=Rc_tetra(:,idx_permutation(kk,1))-Rc_tetra(:,idx_permutation(kk,4));
            B=Rc_tetra(:,idx_permutation(kk,2))-Rc_tetra(:,idx_permutation(kk,4));
            C=Rc_tetra(:,idx_permutation(kk,3))-Rc_tetra(:,idx_permutation(kk,4));
            rc_g2=rc_g-Rc_tetra(:,idx_permutation(kk,4));
            M3=[A B C]'*[A B C];
            N3=[A B C]'*rc_g2;
            cf3=M3\N3;
            d_AB=abs(cf3(3,:)*dot(cross(A,B),C)/norm(cross(A,B)));% distacne to ABD
            M2=[A B]'*[A B];
            N2=[A B]'*rc_g2;
            cf2=M2\N2;
            d_AB(sum(cf2,1)<0|sum(cf2,1)>1)=inf;% distance to surface C
            d_tetra_loop((jj-1)*4+kk,:)=d_AB;
        end
    end

    d_tetra=min(d_tetra_loop,[],1);
    %         d_particle=min(pdist2(rc_g',Rc'),[],2)'-D/2;% 占用内存大

    d1=dot((reshape(rc_g,3,1,[])-Rc0_tmp),repmat(Ori_tmp,1,1,size(rc_g,2)));
    d2=vecnorm((reshape(rc_g,3,1,[])-Rc0_tmp)-d1.*repmat(Ori_tmp,1,1,size(rc_g,2)));
    d1=max(abs(d1)-h/2,0);d2=max(d2-r,0);
    d_particle=min(permute(sqrt(d1.^2+d2.^2),[2,3,1]),[],1);

    r_sphere=min(d_tetra,d_particle);
    clear d_tetra_loop d_AB A B C rc_g2 M3 N3 M2 N2 cf2 surf_common surf_common2 d1 d2 map
    %%
    %     [~,idx_max]=max(r_sphere);
    %     idx_keep=idx_max;
    %     idx_remain=1:size(rc_g,2);
    %     while min(r_sphere(idx_keep))/dx>1
    %         [~,idx_max]=max(r_sphere(idx_remain));
    %         idx_keep=[idx_keep idx_remain(idx_max)];
    %         dist_s=sqrt(sum((rc_g(:,idx_remain)-rc_g(:,idx_remain(idx_max))).^2,1));
    %         idx_remove1=idx_remain(dist_s<=r_sphere(idx_remain(idx_max))*0.7);
    %         idx_remove2=idx_remain(dist_s<=r_sphere(idx_remain(idx_max))*1.1&...
    %             r_sphere(idx_remain)<r_sphere(idx_remain(idx_max))*0.3);
    %         idx_remove=unique([idx_remove1 idx_remove2 idx_remain(idx_max)]);
    %         idx_remain=idx_remain(~ismember(idx_remain,idx_remove));
    %         %         whos idx_remain idx_remove idx_keep
    %     end
    %     LOGIC_draw=false(1,size(rc_g,2));
    %     LOGIC_draw(idx_keep)=1;
    %     %         disp(length(idx_keep))
    %     %         figure(2);clf;hold on
    %     r_sphere1=r_sphere(LOGIC_draw);
    %     rc_g1=rc_g(:,LOGIC_draw);
    r_sphere1=r_sphere;
    rc_g1=rc_g;
    %%
    r_sphere_cell{ii}=r_sphere1;
    rc_g_cell{ii}=rc_g1;
    range_cell{ii}=range;
    %     show_pore(rc_g_cell{ii},r_sphere_cell{ii});
end
rc_g_cell=rc_g_cell(cellfun('length',r_sphere_cell)>1);
range_cell=range_cell(cellfun('length',r_sphere_cell)>1);
r_sphere_cell=r_sphere_cell(cellfun('length',r_sphere_cell)>1);
%%
save
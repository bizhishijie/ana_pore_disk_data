clear;clc
% Rc=Rc(:,1:100);Ori=Ori(:,1:100);
load('p.mat');p=p';n=length(p);
l=0.1;% be same as create lines
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
dx=h/10;
cylinder_size=[r,h];
fileList=dir('../basic/pack/pack*');
Td_map=[2,3,4;1,3,4;1,2,4;1,2,3];
Li_map=[1,2;1,3;1,4;2,3;2,4;3,4];
for nn=1:length(fileList)
    load(['../basic/pack/' fileList(nn).name '/basic.mat'])
    disp(nn)
    cylinder_sample_point=cylinder_to_point(Rc,Ori,p)';
    [Vv,Cv]=voronoin(cylinder_sample_point);
    Td=delaunayn(cylinder_sample_point);
    Td_fix=floor((Td-1)/n)+1;
    rs_pore=zeros(length(Td),3);
    dist_pore=zeros(length(Td),1);
    Nd_idx=cell(length(Td),1);
    Nd_tmp_id=cell(1,length(Td));
    parfor ii=1:length(Td)
        Nd_tmp_id{ii}=find(sum(ismember(Td,Td(ii,:)),2)==3);
    end% neighbors
    %%
    f=@(t,line_tmp,Rc_tmp,Ori_tmp)(-min(distance_point_cylinder((t*(line_tmp(1:3)-line_tmp(4:6))+line_tmp(4:6)),Rc_tmp,Ori_tmp,cylinder_size)));
    parfor ii=1:length(Td)
        %     Cv_tmp=cell2mat(Cv(Td(ii,:))');
        tetra_tmp=cylinder_sample_point(Td(ii,:),:);% 三维的对应的是四面体
        %     plot(tri_tmp(:,1),tri_tmp(:,2),'r')
        DT = delaunayTriangulation(tetra_tmp);
        rs_pore_list = circumcenter(DT);% 四面体的外心

        % Rc_eff_id=Td_fix(ii,:);
        %     Rc_eff_id=unique(Td_fix(ii,:));
        %     Rc_tmp=Rc(:,Rc_eff_id);
        %     Ori_tmp=Ori(:,Rc_eff_id);
        Rc_tmp=Rc;
        Ori_tmp=Ori;
        Rc_eff_id=Td_fix(ii,:);

        if ~is_point_in_tetra(rs_pore_list,tetra_tmp(1,:),tetra_tmp(2,:),tetra_tmp(3,:),tetra_tmp(4,:))
            % 如果外心不在四面体的内部
            rs_pore_list=zeros(4,3);
            dist_pore_list=zeros(4,1);
            is_in_tri=false(4,1);
            is_in_tetra=false(4,1);
            is_in_diff_disk=false(4,1);
            % 如果不在四面体的内部则需要计算四个面的外心的距离
            for jj=1:size(Td_map,1)
                tri_tmp_cir=tetra_tmp(Td_map(jj,:),:);
                C = circumcenter_3d(tri_tmp_cir(1,:),tri_tmp_cir(2,:),tri_tmp_cir(3,:));

                rs_pore_list(jj,:)=C;
                dist_pore_list(jj)=min(distance_point_cylinder(C',Rc_tmp,Ori_tmp,cylinder_size));
                is_in_tri(jj)=is_point_in_tri(C,tri_tmp_cir(1,:),tri_tmp_cir(2,:),tri_tmp_cir(3,:));
                is_in_tetra(jj)=is_point_in_tetra(C,tetra_tmp(1,:),tetra_tmp(2,:),tetra_tmp(3,:),tetra_tmp(4,:));
                is_in_diff_disk(jj)=length(unique(Td_fix(ii,Td_map(jj,:))))>1;
            end
            dist_pore_list=dist_pore_list(is_in_tri&is_in_diff_disk&is_in_tetra);
            if ~isempty(dist_pore_list)
                [dist_pore(ii),max_idx]=max(dist_pore_list);
                rs_pore(ii,:)=rs_pore_list(max_idx,:);
            else% 如果四个三角形也都是钝角三角形
                rs_pore_list=zeros(6,3);
                dist_pore_list=zeros(6,1);
                for jj=1:size(Li_map,1)
                    line_tmp=reshape(tetra_tmp(Li_map(jj,:),:)',[],1);
                    [min_point_position,dist_pore_list(jj)]=fminbnd(@(t)f(t,line_tmp,Rc_tmp,Ori_tmp),0,1);
                    rs_pore_list(jj,:)=min_point_position*(line_tmp(1:3)-line_tmp(4:6))+line_tmp(4:6);
                end
                [dist_pore(ii),max_idx]=max(-dist_pore_list);
                rs_pore(ii,:)=rs_pore_list(max_idx,:);
            end
        else
            % 在四面体内部只需要取这个点就行
            rs_pore_tmp=rs_pore_list;
            dist_pore(ii)=min(distance_point_cylinder(rs_pore_tmp',Rc_tmp,Ori_tmp,cylinder_size));
            rs_pore(ii,:)=rs_pore_tmp;
        end
        %     figure(1);clf;hold on;axis equal
        %     tetramesh(DT)
        %     plot3(rs_pore(ii,1),rs_pore(ii,2),rs_pore(ii,3),'*')
    end
    %%
    % 计算准pore的邻居
    % Nd_tmp_id=cellfun(@(c)find(sum(ismember(Td,c),2)==3),num2cell(Td,2),'UniformOutput',false);\
    dist_pore_tmp=cellfun(@(c)dist_pore(c),Nd_tmp_id,'UniformOutput',false);
    rs_pore_tmp=cellfun(@(c)rs_pore(c,:),Nd_tmp_id,'UniformOutput',false);
    % sum(cellfun('length',Nd_tmp_id)<4);% 邻居比较少的是靠边的
    %% 删掉一部分邻居
    parfor ii=1:length(Td)
        Nd_tmp_loop=Nd_tmp_id{ii};
        dist_pore_loop=dist_pore_tmp{ii};
        rs_pore_loop=rs_pore_tmp{ii};
        Td_fix_tmp=Td_fix(ii,:);
        particle_num=length(unique(Td_fix_tmp));
        % 找邻居的规则：
        % 有三个共同的点的
        % 验证彼此的中心在彼此的球内部
        if particle_num==4||particle_num==3% 如果是在4个颗粒中间的，邻居也是4个的话判断互相在对方球心内
            dist_points=vecnorm(rs_pore(ii,:)-rs_pore_loop,2,2);
            flag1=dist_points<dist_pore(ii)&dist_points<dist_pore_loop;
            flag2=cellfun(@(c)length(unique(c))~=4,num2cell(Td_fix(Nd_tmp_loop,:),2));
            flag3=all(sort(Td_fix_tmp)==sort(Td_fix(Nd_tmp_loop,:),2),2);% 检查是否完全一致
            Nd_idx{ii}=Nd_tmp_loop(flag1|flag2|flag3);
            %     elseif particle_num==3

        elseif particle_num==2
            if sum((Td_fix_tmp-Td_fix_tmp(1))==0)==2% 2+2的模式
                Nd_idx{ii}=Nd_tmp_loop;
            else
                flag1=sum(Td_fix(Nd_tmp_loop,:)-mode(Td_fix_tmp)==0,2)~=3;
                flag2=all(ismember(Td_fix(Nd_tmp_loop,:),Td_fix_tmp),2);
                Nd_idx{ii}=Nd_tmp_loop(flag1|flag2);
            end
        else % 只由一个颗粒所定义，只在部分堆积外部中具有，计算精度问题
            %         disp(ii)
            Nd_idx{ii}=Nd_tmp_loop;
            %         Nd_idx{ii}=[];
        end
        %     if isempty(Nd_idx{ii})
        %         disp(1)
        %     end
    end
    % plot(cellfun('length',Nd_idx))
    % axis([0 1 0 1])
    % figure(2)
    % clf
    % hold on
    % edge_pore=[min(line_sample_point,[],1);max(line_sample_point,[],1)];
    % edge_rate=0.1;
    % edge_pore=[edge_pore(:,1)+(edge_pore(:,2)-edge_pore(:,1))*edge_rate edge_pore(:,2)-(edge_pore(:,2)-edge_pore(:,1))*edge_rate];
    % pore_inside=find(all(rs_pore>edge_pore(1,:)&rs_pore
    %             plot(line_pore(:,[1,3]),line_pore(:,[2,4]),'g')<edge_pore(2,:),2));
    %% 第一次合并
    IDX_merge=cell(length(Td),1);
    dist_pore_tmp2=cellfun(@(c)dist_pore(c),Nd_idx,'UniformOutput',false);
    Td_pore_tmp2=cellfun(@(c)dist_pore(c),Nd_idx,'UniformOutput',false);
    %%
    parfor ii=1:size(Td,1)
        n_idx0=Nd_idx{ii};
        dist_tmp=dist_pore_tmp2{ii};
        n_idx1=n_idx0(dist_tmp==max(dist_tmp)&dist_tmp>dist_pore(ii));
        % if length(n_idx1)>1
        %     n_idx1=n_idx1(1);
        % end
        n_idx2=n_idx0(length(unique(Td_fix(ii,:)))>2&all(sort(Td_fix(ii,:))==sort(Td_fix(n_idx0,:),2),2));
        IDX_merge{ii}=unique([n_idx1;n_idx2]);
    end% 每个三角形合并的三角形的id
    IDX_merge1 = cell_combine_(IDX_merge);
    % plot(cellfun('length',IDX_merge1))
    %% 第二次合并
    idx_tmp=1:length(IDX_merge1);
    idx_tmp=repelem(idx_tmp,cellfun('length',IDX_merge1));
    IDX_convert(cell2mat(IDX_merge1))=idx_tmp;% 每个四面体中心merge到的编号
    IDX_merge_pore=cell(1,length(IDX_merge1));
    merge_particle_cell=cell(1,length(IDX_merge1));
    parfor ii=1:length(IDX_merge1)% 第二次合并
        IDX_merge1_tmp=IDX_merge1{ii};
        merge_particle=unique(Td_fix(IDX_merge1_tmp,:));% 这个pore对应的颗粒
        merge_particle_cell{ii}=merge_particle;
        pore_neighbor=setdiff(unique(cell2mat(Nd_idx(IDX_merge1_tmp))),IDX_merge1_tmp);% pore的邻居
        merge_neighbor=IDX_convert(pore_neighbor);% merge的邻居
        %     cellfun(@(c)unique(Td_fix(IDX_merge1{c})),num2cell(merge_neighbor),'UniformOutput',false);
        %     merge_neighbor_particle=cellfun(@(c)unique(c),num2cell(Td_fix(merge_neighbor,:),2),'UniformOutput',false);
        merge_neighbor_particle=cellfun(@(c)unique(Td_fix(c,:)),IDX_merge1(merge_neighbor),'UniformOutput',false);
        is_neighbor=cellfun(@(c)all(ismember(merge_particle,c)),merge_neighbor_particle);
        %     any(is_neighbor),length(merge_particle)
        pore_neighbor=pore_neighbor(is_neighbor);
        merge_neighbor=merge_neighbor(is_neighbor);
        [~,max_idx]=max(dist_pore(pore_neighbor));
        IDX_merge_pore{ii}=merge_neighbor(max_idx);
    end

    IDX_merge_pore_1 = cell_combine_(IDX_merge_pore);
    IDX_merge2=cellfun(@(c)cell2mat(IDX_merge1(c)),IDX_merge_pore_1,'UniformOutput',false);
    %%
    % merge_particle_cell2=zeros(1,length(IDX_merge2));
    % parfor ii=1:length(IDX_merge2)
    %     IDX_merge2_tmp=IDX_merge2{ii};
    %     merge_particle=unique(Td_fix(IDX_merge2_tmp,:));% 这个pore对应的颗粒
    %     merge_particle_cell2(ii)=length(merge_particle);
    % end

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
        Ori_tmp=Ori(:,unique(Td_fix(idx_tetra,:)));
        %         if any(~is_in(idx_dcell))
        %             continue
        %         end
        Rc_d=cylinder_sample_point(idx_dcell,:)';
        range=[min(Rc_d,[],2) max(Rc_d,[],2)];
        [xg,yg,zg]=meshgrid(range(1,1):dx:range(1,2),range(2,1):dx:range(2,2),range(3,1):dx:range(3,2));
        rc_g=[xg(:) yg(:) zg(:)]';
        if length(rc_g)>1000000% 超出内存限制，多半是没有合并好
            continue
        end
        clear Rc_d idx_dcell xg yg zg
        %% 判断点是否需要留下
        d1=dot((reshape(rc_g,3,1,[])-Rc_tmp),repmat(Ori_tmp,1,1,size(rc_g,2)));
        d2=vecnorm((reshape(rc_g,3,1,[])-Rc_tmp)-d1.*repmat(Ori_tmp,1,1,size(rc_g,2)));
        d1=max(abs(d1)-h*1.1/2,0);d2=max(d2-r,0);% d1需要减去大一些以填补缝隙
        idx_out_particle=permute(all(d1>0|d2>0,[1,2]),[2,3,1]);% 在同一个圆盘的两个面的外面
        clear d1 d2
        idx_in_tetra_loop=false(length(idx_tetra),size(rc_g,2));
        %%
        parfor jj=1:length(idx_tetra)
            Rc_tetra=cylinder_sample_point(Td(idx_tetra(jj),:),:)';
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
        idx_dcell=unique(Td_fix(idx_tetra,:));

        rc_g_cell{ii}=rc_g;
        idx_dcell_cell{ii}=idx_dcell;
    end
    load(['../basic/pack/' fileList(nn).name '/edge.mat'])
    range=[max(edge,[],2),min(edge,[],2)];
    cell_not_empty=~cellfun('isempty',rc_g_cell);
    rc_g_cell=rc_g_cell(cell_not_empty);
    idx_dcell_cell=idx_dcell_cell(cell_not_empty);
    cell_in=cellfun(@(c)all(~isempty(c)&c<range(:,1)&c>range(:,2),"all"),rc_g_cell);
    cell_eff=cellfun('length',idx_dcell_cell)>3;
    rc_g_cell=rc_g_cell(cell_in&cell_eff);
    idx_dcell_cell=idx_dcell_cell(cell_in&cell_eff);
    % clearvars -except rc_g_cell idx_dcell_cell nn fileList
    save(['../basic/pack/' fileList(nn).name '/pore_data.mat'],'rc_g_cell','idx_dcell_cell')
end
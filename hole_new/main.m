% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
load('p.mat');p=p';p_length=length(p);
for nn=1:length(fileList)
    %% 加载数据并计算voronoi
    load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
    load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
    load(['..\basic\pack\' fileList(nn).name '\edge.mat'])
    disp(fileList(nn).name)
    n=length(Rc);
    %     n=100;
    p_all=cylinder_to_point(Rc,Ori,p);
    p_all=p_all(:,1:p_length*n);
    [v,c]=voronoin(p_all');
    %% 计算每个颗粒的采样点对应的连接
    v_cell=cell(1,n);% 保存每个颗粒的元胞的边界点
    k_cell=cell(1,n);% 保存每个颗粒的元胞的边界点的连接
    parfor pp=1:n% 改成全部粒子
        ctmp=c((pp-1)*p_length+1:p_length*pp);
        vp=cell2mat(ctmp');
        tbl=tabulate(vp);
        vb=tbl(tbl(:,2)==1|tbl(:,2)==2|tbl(:,2)==3,1);
        v_cell{pp}=vb;
        % 第pp个颗粒的 voronoi 的外侧边界的点的编号
        k_t=cell(1,p_length);
        for ii=1:p_length
            p_id=ctmp{ii};
            P=v(p_id,:);
            if ~any(any(isinf(P)))
                k=boundary(P,0);
                k_t{ii}=p_id(k);% 转化成v的编号
            else
                k=boundary(P(2:end,:),0);% 去掉无穷大的
                k_t{ii}=p_id(k+1);% 转化成v的编号
            end

            %             plot3(P(:,1),P(:,2),P(:,3),'.','MarkerSize',10)
            %             hold on
            %             trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor','red','FaceAlpha',0.1)

        end
        k_cell{pp}=k_t;% 保存了每个颗粒的采样点的边界的连接方式
    end
    %% 计算每个圆盘相邻的圆盘
    ne=cell(1,n);
    for pp=1:n
        ne{pp}=[];
    end
    for pp=1:n
        for qq=pp+1:n
            cp=intersect(v_cell{pp},v_cell{qq});
            if ~isempty(cp)
                ne{pp}(end+1)=qq;
                ne{qq}(end+1)=pp;
            end
        end
    end
    %% 计算孔的位置和孔的rk属性
    tbl=tabulate(cell2mat(v_cell'));
    map=find(tbl(:,2)==4);
    pore_rc=v(map,:);
    pore_is_in=false(1,length(pore_rc));
    co_edge=convhulln(edge');
    for ii=1:length(pore_rc)
        pore_is_in(ii)=inpolygon3d(pore_rc(ii,:)',edge,co_edge);
    end
    % 删掉不在 edge 里面的孔
    pore_rc=pore_rc(pore_is_in,:);
    map=map(pore_is_in);
    rk=zeros(1,length(pore_rc));
    for ii=1:length(pore_rc)
        dis_tmp=distance_point_cylinder(pore_rc(ii,:)',Rc,Ori,cylinder_size);
        rk(ii)=min(dis_tmp(1,:));
    end
    pore_rc=pore_rc(rk>0,:);
    map=map(rk>0);% 孔对应的点在v中的编号
    rk=rk(rk>0);
    c_pore_id=cellfun(@(v1)map(ismember(map,v1(v1~=1))),v_cell,'UniformOutput',false);
    % 去除掉第一个的inf
    % 每个圆盘附近的pore的id
    %% 计算每个采样点相邻的圆盘
    co=cell(1,n);% 表示颜色
    parfor pp=1:n
        co{pp}=cell(1,length(k_cell{pp}));
        for jj=1:length(k_cell{pp})
            disp(jj)
            vk=k_cell{pp}{jj};% 这个粒子上的第jj个采样点
            co{pp}{jj}=zeros(1,size(vk,1));
            io=zeros(1,size(vk,1));% is outside
            for kk=1:size(vk,1)
                v_tmp=vk(kk,:);% 该采样点的第kk个三角面
                nc=cell2mat(cellfun(@(v)any(all(ismember(v_tmp,v))),v_cell,'UniformOutput',false));
                % 如果这三个采样点全都属于某一个粒子的表面，则标记
                nc(pp)=0;% 去除掉自己
                if any(nc)
                    co{pp}{jj}(kk)=find(nc,1);
                end
            end
        end
    end
    %% 取消注释画颗粒和元胞
    pp=30;
    show_particle(Rc(:,pp),Ori(:,pp),v,k_cell{pp},v(c_pore_id{pp},:),co{pp})
    l=light;
    l.Position=[1,-1,1];
    axis off
    

    %% 计算pore的诸多参数
    pore_c_id=zeros(size(pore_rc,1),4);% 每个孔相邻的圆柱
    for pp=1:size(pore_rc,1)
        pore_c_id(pp,:)=find(cellfun(@(v)ismember(map(pp),v),v_cell));
    end

    neighbor_cell=cell(1,size(pore_rc,1));% 孔的邻居
    for pp=1:size(pore_rc,1)
        neighbor_cell{pp}=find(sum(ismember(pore_c_id,pore_c_id(pp,:)),2)==3);
    end
    %%
    path_cell=cell(1,size(pore_rc,1));% 孔的路径
    throat_cell=cell(1,size(pore_rc,1));% 嗓子的位置
    parfor pp=1:size(pore_rc,1)
        path_cell{pp}=cell(1,length(neighbor_cell{pp}));
        throat_cell{pp}=cell(1,length(neighbor_cell{pp}));
        for jj=1:length(neighbor_cell{pp})% 对每个孔和它的邻居遍历
            qq=neighbor_cell{pp}(jj);
            cc=intersect(pore_c_id(pp,:),pore_c_id(qq,:));% common cylinder
            tbl=tabulate(cell2mat(v_cell(cc)'));
            path_tmp=tbl(tbl(:,2)==3,1);% 需要根据远近将path排序，便于合并时连接

            dis=sqrt(sum((v(path_tmp,:)-pore_rc(pp,:)).^2,2));
            path_tmp=sortrows([dis,path_tmp]);
            path_tmp=path_tmp(:,2);
            path_cell{pp}{jj}=path_tmp;% 排序后放入cell

            Rc_t=Rc(:,cc);Ori_t=Ori(:,cc);
            dl=zeros(1,length(path_tmp));% 根据距离最短找到嗓子的位置
            for kk=1:length(path_tmp)
                pa=path_tmp(kk);
                dis=distance_point_cylinder(v(pa,:),Rc_t,Ori_t,cylinder_size);
                dis=min(dis(1,:));
                dl(kk)=dis;
            end
            [~,idx]=min(dl);
            dp=diff(v(path_tmp,:));
            if idx>size(dp,1)% 如果第一个点就是throat，实际上可能是由于采样点数目不够多造成的
                dp=dp(idx-1,:);
            else
                dp=dp(idx,:);
            end
            throat_cell{pp}{jj}=[path_tmp(idx) dp/norm(dp)];
            % 第一个数是嗓子所在点的编号，后三个数是path在嗓子处的朝向，朝向总是向外的
        end
    end
    %% 计算孔拥有的采样点
    ks=cell(1,size(pore_rc,1));% 孔有的采样点和连接
    f=@(p1,p0,ori)sign(sum(repmat(ori,size(p1,1),1).*(p1-p0),2));% 判断一个点在平面的哪个方向的函数
    parfor pp=1:size(pore_rc,1)
        ks{pp}=[];
        for jj=1:4
            ctmp=pore_c_id(pp,jj);% cylinder
            stmp= k_cell{ctmp}(ismember(co{ctmp},pore_c_id(pp,:)));% 孔在这个圆盘上拥有的采样点
            for kk=1:length(stmp)
                p_id=unique(stmp{kk}(:));
                sutmp=v(p_id,:);% 获得了边界点，接下来需要根据throat的位置和方向删除
                flag=true(size(sutmp,1),1);
                for tt=1:length(throat_cell{pp})
                    s_2=f(sutmp,v(throat_cell{pp}{tt}(1),:),throat_cell{pp}{tt}(2:4));
                    s_1=f(pore_rc(pp,:),v(throat_cell{pp}{tt}(1),:),throat_cell{pp}{tt}(2:4));
                    flag=flag & (s_1~=s_2);
                end
                p_id=p_id(flag,:);% 找到了在内侧的点的id
                s=stmp{kk}(all(ismember(stmp{kk},p_id),2),:);
                if ~isempty(s)
                    stmp{kk}=s;
                end
            end
            ks{pp}=[ks{pp} stmp];
        end
    end
    %% 展示pore
    pp=20;
    show_pore(pore_rc(pp,:),ks{pp},v,Rc(:,neighbor_cell{pp}),Ori(:,neighbor_cell{pp}),pore_rc(neighbor_cell{pp},:),path_cell(pp));
end
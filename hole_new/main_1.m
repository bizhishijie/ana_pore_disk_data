% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
load('p.mat');p=p';p_length=length(p);
for nn=1:15
    disp(nn)
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
    parfor pp=1:n
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
            if ~any(isinf(P),'all')
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
    map=tbl(tbl(:,2)==4,1);
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
    %     pore_rc=pore_rc(rk>0,:);
    %     map=map(rk>0);% 孔对应的点在v中的编号
    %     rk=rk(rk>0);
    c_pore_id=cellfun(@(v1)map(ismember(map,v1(v1~=1))),v_cell,'UniformOutput',false);
    % 去除掉第一个的inf
    % 每个圆盘附近的pore的id
    %% 计算每个采样点相邻的圆盘
    co=cell(1,n);% 表示颜色
    parfor pp=1:n
        %         disp(pp)
        co{pp}=cell(1,length(k_cell{pp}));
        for jj=1:length(k_cell{pp})
            vk=k_cell{pp}{jj};% 这个粒子上的第jj个采样点
            co{pp}{jj}=zeros(1,size(vk,1));
            io=zeros(1,size(vk,1));% is outside
            for kk=1:size(vk,1)
                v_tmp=vk(kk,:);% 该采样点的第kk个三角面
                nc=ne{pp}(cell2mat(cellfun(@(v)any(all(ismember(v_tmp,v))),v_cell(ne{pp}),'UniformOutput',false)));
                % 如果这三个采样点全都属于某一个粒子的表面，则标记
                if any(nc)
                    co{pp}{jj}(kk)=nc(1);
                end
            end
        end
    end
    save(['..\basic\pack\' fileList(nn).name '\voro.mat'],'v','pore_rc','co','k_cell','v_cell','map','rk');
end
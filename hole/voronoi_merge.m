% 合并一个圆盘上的采样点的格子
clear
d=30.75*2;
h=4.81*2;
r=d/2;
cylinder_size=[r,h];
fileList=dir('..\basic\pack\pack*');
load('tri.mat')
p_length=size(p,2);
parfor ii=7:length(fileList)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Rc').Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat'],'Ori').Ori;
    disp(fileList(ii).name)

    regions=load(['..\basic\pack\' fileList(ii).name '\voronoi.mat'],'regions').regions;
    vertice=load(['..\basic\pack\' fileList(ii).name '\voronoi.mat'],'vertice').vertice;

    dt_cell=cell(1,size(Rc,2));
    volume_list=zeros(1,size(Rc,2));

    Rc_min=min(Rc,[],2)+r;
    Rc_max=max(Rc,[],2)-r;
    box=[Rc_min Rc_max];

    for jj=1:size(Rc,2)
        dt_tmp=[];
        volume=0;
        % 对每一个圆盘遍历，对应regions中的前若干个点
        id_tmp=regions((jj-1)*p_length+1:jj*p_length);
        id_tmp=cell2mat(cellfun(@transpose,id_tmp,'UniformOutput',false));% 获得全部的id列表
        id_tmp=sort(id_tmp);
        id_uni=unique(id_tmp);
        id_statistics=zeros(length(id_uni),2);
        id_statistics(:,1)=id_uni;% 统计每一个点出现的频次
        for kk=1:size(id_uni)
            id_statistics(kk,2)=sum(id_tmp==id_uni(kk));
        end
        id_uni(id_statistics(:,2)>=4)=[];
        % 在 id  _tmp 里面，出现4次及以上的是内部的点，应该删除
        cylinder_point=vertice(id_uni,:);
        % 找到对应的点的编号
        if ~any(any(isinf(cylinder_point)))% 如果圆盘不在最外侧
            if is_in_box(cylinder_point',box)% 如果圆盘的格子不超过限制
                [~,volume] = convhull(cylinder_point,"Simplify",true);
                dt_tmp=delaunayTriangulation(cylinder_point);
            end
        else
            continue
        end
        %     defaultFaceColor  = [0.6875 0.8750 0.8984];
        %     trisurf(edge, cylinder_point(:,1),cylinder_point(:,2),cylinder_point(:,3) , ...
        %         'FaceColor', defaultFaceColor, 'FaceAlpha',0.8)
        %     hold on
        %     disp(jj)
        dt_cell{jj}=dt_tmp;
        volume_list(jj)=volume;
    end
    parsave(['..\basic\pack\' fileList(ii).name '\voronoi_cut_dt.mat'],dt_cell);
    parsave(['..\basic\pack\' fileList(ii).name '\voronoi_volume.mat'],volume_list);
end
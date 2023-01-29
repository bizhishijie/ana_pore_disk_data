% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];%
s=20*d;
rc_num=1e5;
parfor ff=1:length(fileList)
    warning('off')
    Rc=load(['..\basic\pack\' fileList(ff).name '\basic.mat']).Rc;
    del_list=[];
    for ii=1:length(Rc)
        dis=sqrt( sum((Rc-Rc(:,ii)).^2) );
        if min( dis(dis~=0) )>1.3*r
            del_list=[del_list ii];
        end
    end
    % 删除离其他圆盘过远的孤立圆盘
    Rc0=Rc;
    Rc(:,del_list)=[];
    disp(length(del_list));

    disp(fileList(ff).name)
    K=convhull(Rc');
    Rc_mean=mean(Rc,2);
    rc=[0 1 1;0 1 -1; 0 -1 1; 0 -1 -1;
        2 1 1;2 1 -1; 2 -1 1; 2 -1 -1];

    edge=(rc)'*20*d;
    clf
    for ii = 1:length(K)
        p1=Rc(:,K(ii,1));
        p2=Rc(:,K(ii,2));
        p3=Rc(:,K(ii,3));

        v1=p1-p2;
        v2=p1-p3;
        v1=v1/norm(v1);
        v2=v2/norm(v2);
        n=cross(v1,v2);
        n=n/norm(n)*sign(dot(n,Rc_mean-p1));
        % 修正法向量的方向，修正后方向是向内的
        flate_v=null(n');
        vf1=flate_v(:,1);
        vf2=flate_v(:,2);
        % 根据零空间获得平面的两个向量
        rc_tmp=zeros(size(rc));
        base=[n vf1 vf2]';% 转置后每一行是一个基向量
        for jj =1:length(rc)
            rc_tmp(jj,:)=rc(jj,:)*base;
        end
        %     rc_tmp=rc_tmp./sqrt(sum(rc_tmp.^2,2));
        [~,edge1]=intersects_polyhedron_(edge,(rc_tmp'*s)+p1+n*r);
        % 新生成的和原来的求交集
        if ~isempty(edge1)
            edge=edge1;
        end
    end
    parsave(['..\basic\pack\' fileList(ff).name '\edge.mat'],edge)
    k=convhulln(edge');
    h=trisurf(k,edge(1,:),edge(2,:),edge(3,:),ones(1,size(edge,2)),'faceAlpha',0.5);
    hold on
    %     show_cylinder(Rc,Ori)
    plot3(Rc(1,:),Rc(2,:),Rc(3,:),'o')
    axis equal
    saveas(gcf,['..\basic\pack\' fileList(ff).name '\edge.jpg'])
    saveas(gcf,['..\basic\pack\' fileList(ff).name '\edge.fig'])
    % 保存每个圆盘的质心是否在中间的区域
    is_in=false(1,length(Rc0));
    for jj=1:length(Rc0)
        if ~ismember(jj,del_list) && inpolygon3d_new(Rc0(:,jj),edge)
            is_in(jj)=true;
        end
    end

    parsave(['..\basic\pack\' fileList(ff).name '\rc_is_in.mat'],is_in)
end
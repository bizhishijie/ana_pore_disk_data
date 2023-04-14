% 读取圆盘位置
clear
fileList=dir('../basic/data/*_plate');
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%
rc_num=1e6;
for ii=1:length(fileList)
    rc_list=zeros(2,rc_num);
    load(['../basic/data/' fileList(ii).name '/basic_adjust.mat'])
    load(['../basic/data/' fileList(ii).name '/edge.mat'])

    idx_keep=setxor(1:size(Rc,2),idx_remove);
    Rc=[Rc(:,idx_keep) Rc_add];
    Ori=[Ori(:,idx_keep) Ori_add];

    disp(fileList(ii).name)
    rand_point=zeros(3,rc_num);

    %     Rc_max=max(Rc,[],2)-r;
    %     Rc_min=min(Rc,[],2)+r;

    Rc_max_sphere=max(edge,[],2);
    Rc_min_sphere=min(edge,[],2);
    %
    %     %     Rc_outside=any((Rc>Rc_max)|(Rc<Rc_min));
    %     Rc_1=Rc(:,inside_id);
    %     Ori_1=Ori(:,inside_id);% 去除比较靠外的圆盘
    %
    %     Rc_2=Rc(:,all(Rc>Rc_min_sphere&Rc<Rc_max_sphere));
    %     dt=delaunayTriangulation(Rc_2');
    %     poly=convexHull(dt);% 质心组成的多边形
    %     polygon=Rc_2(:,unique(poly(:)));%

    % 取的小球不太靠边
    parfor sphere_id=1:rc_num
        sphere_Rc=rand(3,1).*(Rc_max_sphere-Rc_min_sphere)+Rc_min_sphere;
        while is_point_in_cylinder(sphere_Rc,Rc,Ori,cylinder_size)||(~inpolygon3d_new(sphere_Rc,edge))
            % 如果随机到的点在圆盘里面或者随机到的点在多面体外面
            sphere_Rc=rand(3,1).*(Rc_max_sphere-Rc_min_sphere)+Rc_min_sphere;
        end
        rc_tmp_list=[distance_point_cylinder(sphere_Rc,Rc,Ori,cylinder_size);Rc;Ori];
        rc_tmp_list=sortrows(rc_tmp_list')';
        rc_list(:,sphere_id)=rc_tmp_list(1:2,1);% 取最近的接触点的
        rand_point(:,sphere_id)=sphere_Rc;
        %         Rc_tmp=rc_tmp_list(3:5,1:20);
        %         Ori_tmp=rc_tmp_list(6:8,1:20);
        %         show_cylinder(Rc_tmp,Ori_tmp);
        %         [X,Y,Z] =sphere(20);
        %         X=X*rc_tmp_list(1)+sphere_Rc(1);
        %         Y=Y*rc_tmp_list(1)+sphere_Rc(2);
        %         Z=Z*rc_tmp_list(1)+sphere_Rc(3);
        %         surf(X,Y,Z)
        %         axis equal

    end
    save(['../basic/data/' fileList(ii).name '/rc_list.mat'],'rc_list')
    %     save(['../basic/pack/' fileList(ii).name '/rand_point.mat'],'rand_point')
end

% histogram(rc_list(1,:));
% hold on
% rc_list_1=rc_list(:,rc_list(2,:)==1);
% rc_list_2=rc_list(:,rc_list(2,:)==2);
% rc_list_3=rc_list(:,rc_list(2,:)==3);
% [counts_1,centers_1] = hist(rc_list_1(1,:),50);
% [counts_2,centers_2] = hist(rc_list_2(1,:),50);
% [counts_3,centers_3] = hist(rc_list_3(1,:),50);
% [counts_a,centers_a] = hist(rc_list(1,:),50);
% plot(centers_1,counts_1)
% plot(centers_2,counts_2)
% plot(centers_3,counts_3)
% plot(centers_a,counts_a)
% set(gca,'yScale','log')
% legend('角','底','面','总')
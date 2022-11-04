% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
d=30.75*2;
h=4.81*2;
r=d/2;
cylinder_size=[r,h];

for ii=1:length(fileList)
    rc_list=[];
    load(['..\basic\pack\' fileList(ii).name '\basic.mat'])
    disp(fileList(ii).name)

    Rc_max=max(Rc,[],2)-r;
    Rc_min=min(Rc,[],2)+r;

    Rc_outside=find(any((Rc>Rc_max)|(Rc<Rc_min)));
    Rc(:,Rc_outside)=[];
    Ori(:,Rc_outside)=[];% 去除比较靠外的圆盘
    % 取的小球不太靠边
    for sphere_id=1:1e4
        sphere_Rc=rand(3,1).*(Rc_max-Rc_min)+Rc_min;
        while is_point_in_cylinder(sphere_Rc,Rc,Ori,cylinder_size)
            sphere_Rc=rand(3,1).*(Rc_max-Rc_min)+Rc_min;
        end
        rc_tmp_list=[distance_point_cylinder(sphere_Rc,Rc,Ori,cylinder_size);Rc;Ori];
        rc_tmp_list=sortrows(rc_tmp_list')';
        rc_list=[rc_list rc_tmp_list(1:2,1)];% 取最近的接触点的

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
    save(['..\basic\pack\' fileList(ii).name '\rc_list.mat'],'rc_list')
end
histogram(rc_list(1,:));
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
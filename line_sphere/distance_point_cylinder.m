function dis=distance_point_cylinder(point,cylinder_Rc,cylinder_Ori,cylinder_size)
% point 是点的坐标
% cylinder_Rc 是圆柱质心的坐标
% cylinder_Ori 是圆柱顶/底的法向量的朝向
% cylinder_size 是圆柱的尺寸 [r,h]

% d=30.4*2;
% h=4.7*2;
% r=d/2;
% load('..\basic\pack\pack104\basic.mat')
% distance_point_cylinder([1;1;1],Rc(:,1),Ori(:,1),[r,h])

r=cylinder_size(1);
h=cylinder_size(2);
dis=zeros(2,size(cylinder_Rc,2));

RP=point-cylinder_Rc;% 圆盘质心到点的矢量
RP_z_length=dot(RP,cylinder_Ori);% 轴向分量
RP_z=RP_z_length.*cylinder_Ori;
RP_r=RP-RP_z;
RP_r_length=sqrt(sum(RP_r.^2));
RP_z_length=abs(RP_z_length);
% 分三种情况讨论
case_1=RP_z_length>h/2&RP_r_length>r;
case_2=RP_z_length>h/2&RP_r_length<r;
case_3=RP_z_length<h/2&RP_r_length>r;

dis(:,case_1)=[sqrt((RP_z_length(case_1)-h/2).^2+(RP_r_length(case_1)-r).^2);ones(1,sum(case_1))*1];
dis(:,case_2)=[RP_z_length(case_2)-h/2;ones(1,sum(case_2))*2];
dis(:,case_3)=[RP_r_length(case_3)-r;ones(1,sum(case_3))*3];

% dis=(sqrt((RP_z_length-h/2).^2+(RP_r_length-r).^2)).*(RP_z_length>h/2&RP_r_length>r)+...
%     (RP_z_length-h/2).*(RP_z_length>h/2&RP_r_length<r)+...
%     (RP_r_length-r).*(RP_z_length<h/2&RP_r_length>r);
end
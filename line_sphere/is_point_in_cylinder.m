function is_in=is_point_in_cylinder(point,cylinder_Rc,cylinder_Ori,cylinder_size)
% 判断点是否在圆盘内部
% point 是点的坐标
% cylinder_Rc 是圆柱质心的坐标
% cylinder_Ori 是圆柱顶/底的法向量的朝向
% cylinder_size 是圆柱的尺寸 [r,h]
% true 表示点在圆盘内部
r=cylinder_size(1);
h=cylinder_size(2);
dis_1=distance_point_line(point,cylinder_Rc,cylinder_Ori);
dis_2=sqrt(sum((point-cylinder_Rc).^2)-dis_1.^2);
is_in=any(dis_1<r&dis_2<h/2);
end 
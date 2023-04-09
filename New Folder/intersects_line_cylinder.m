function cross_p=intersects_line_cylinder(line_p,line_forward,cylinder_Rc,cylinder_Ori,cylinder_size)
% 返回直线和圆柱相交的线段的交点
r=cylinder_size(1);
h=cylinder_size(2);

% 直线与圆柱面的交点
cross_p_1=intersects_line_cylinder_surface(line_p,line_forward,cylinder_Rc,cylinder_Ori,cylinder_size);
cross_p_1_1=cross_p_1(:,1);
cross_p_1_2=cross_p_1(:,2);
flag_1_1=sqrt(norm(cross_p_1_1-cylinder_Rc)^2-r^2)>h/2||any(isnan(cross_p_1_1));
flag_1_2=sqrt(norm(cross_p_1_2-cylinder_Rc)^2-r^2)>h/2||any(isnan(cross_p_1_1));
% 如果交点超出了圆柱的范围认为无效

% 直线与上下平面的交点
cross_p_2_1=intersects_line_flat(line_p,line_forward,cylinder_Rc+cylinder_Ori*h/2,cylinder_Ori);
cross_p_2_2=intersects_line_flat(line_p,line_forward,cylinder_Rc-cylinder_Ori*h/2,cylinder_Ori);
flag_2_1=sqrt(norm(cross_p_2_1-cylinder_Rc)^2-(h/2)^2)>r;
flag_2_2=sqrt(norm(cross_p_2_2-cylinder_Rc)^2-(h/2)^2)>r;
% 如果交点超出了圆柱的范围认为无效

% 取四个交点中的两个，此距离即为直线和圆柱相交的线段的长度
cross_p=[];
if ~flag_1_1
    cross_p=[cross_p cross_p_1_1];
end
if ~flag_1_2
    cross_p=[cross_p cross_p_1_2];
end
if ~flag_2_1
    cross_p=[cross_p cross_p_2_1];
end
if ~flag_2_2
    cross_p=[cross_p cross_p_2_2];
end
end
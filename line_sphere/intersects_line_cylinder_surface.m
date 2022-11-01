function cross_p=intersects_line_cylinder_surface(line_p,line_forward,cylinder_Rc,cylinder_Ori,cylinder_size)
% 返回直线和圆柱面的交点
r=cylinder_size(1);
h=cylinder_size(2);
if all(line_forward~=cylinder_Ori)||all(line_forward~=-cylinder_Ori)
    [p_t_line,p_t_cylinder]=distance_line_line(line_p,line_forward,cylinder_Rc,cylinder_Ori);
    dist=norm(p_t_line-p_t_cylinder);% 直线到圆柱中轴线的最近距离
    if r<dist
        cross_p=nan(3,2);
        % 和圆没有交点
    else
        chord_length_half=sqrt(r^2-dist^2);% 半弦长
        cos_th=dot(line_forward,cylinder_Ori);
        sin_th=sqrt(1-cos_th^2);
        v_cross_dist=chord_length_half/sin_th*line_forward;% 距离最近处到交点的向量
        cross_p=[p_t_line+v_cross_dist p_t_line-v_cross_dist];
    end
else
    cross_p=nan(3,2);
end
end
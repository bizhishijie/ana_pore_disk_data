function [p_t1,p_t2]=distance_line_line(line_1_p,line_1_forward,line_2_p,line_2_forward)
% 计算异面曲线的最近点
line_1_forward = line_1_forward/norm(line_1_forward);
line_2_forward = line_2_forward/norm(line_2_forward);
flat_normal = cross(line_1_forward,line_2_forward);
p_t1 = (cross(flat_normal,line_2_forward)'*(line_2_p-line_1_p))/(cross(flat_normal,line_2_forward)'*line_1_forward)*line_1_forward + line_1_p;
p_t2 = (cross(flat_normal,line_1_forward)'*(line_1_p-line_2_p))/(cross(flat_normal,line_1_forward)'*line_2_forward)*line_2_forward + line_2_p;
end
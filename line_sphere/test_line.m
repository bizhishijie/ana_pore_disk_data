% 检测圆盘和直线的交点的计算是否正确
% 读取圆盘位置
clf
load('..\basic\pack\pack6\basic.mat')

d=30.75*2;
h=4.81*2;
r=d/2;
cylinder_size=[r,h];

Rc_tmp=Rc(:,1);
Ori_tmp=Ori(:,1);
Rc_max=max(Rc,[],2);
Rc_min=min(Rc,[],2);
% 生成直线
% line_p=rand(3,1).*(Rc_max-Rc_min)+Rc_min;
% line_forward=rand(3,1);
% line_forward=line_forward/norm(line_forward);
line_p=[460;462;430];
line_forward=[1;-1;0.8];
line_forward=line_forward/norm(line_forward);

n_length=200;% 画的直线的长度
draw_cylinder(Rc_tmp,Ori_tmp,'r',r,h,1)
p_1=line_p-n_length*line_forward;
p_2=line_p+n_length*line_forward;
line([p_1(1) p_2(1)],[p_1(2) p_2(2)],[p_1(3) p_2(3)])
axis equal

cross_p=intersects_line_cylinder(line_p,line_forward,Rc_tmp,Ori_tmp,cylinder_size);
for ii=1:size(cross_p,2)
    plot3(cross_p(1,ii),cross_p(2,ii),cross_p(3,ii),'b*');
end
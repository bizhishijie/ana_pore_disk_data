load('..\basic\pack\pack104\basic.mat')
load('G:\毕个业\basic\pack\pack104\rc_is_in.mat')
load('..\basic\pack\pack104\point_pore.mat')
point_pore=point_pore';
Rc=Rc(:,is_in);
Ori=Ori(:,is_in);
d=30.4*2/0.8;
h=4.7*2/0.8;
r=d/2;
cylinder_size=[r,h];

% plot3(Rc(1,:),Rc(2,:),Rc(3,:),'o')
% hold on
% plot3(point_pore(1,:),point_pore(2,:),point_pore(3,:),'o');
% axis equal
% 对于孔的点，到周边的4个圆柱的距离相等
for ii=1:length(point_pore)
    point=point_pore(:,ii);
    dis=distance_point_cylinder(point,Rc,Ori,cylinder_size);
    dis=dis(1,:);
    near_id=find(dis==min(dis));
end
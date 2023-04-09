function show_cylinder(Rc,Ori,c)
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
for ii=1:size(Rc,2)
    draw_cylinder(Rc(:,ii),Ori(:,ii),c,r,h,0.5)
end
axis equal
end
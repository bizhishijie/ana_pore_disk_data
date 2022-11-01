function show_cylinder(Rc,Ori)
d=30.75*2;
h=4.81*2;
r=d/2;
for ii=1:size(Rc,2)
    draw_cylinder(Rc(:,ii),Ori(:,ii),'r',r,h,0.5)
end
axis equal
end
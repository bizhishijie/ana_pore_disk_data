function show_sphere(Rc)
r=23.7/2;
[x,y,z]=sphere();
for ii=1:size(Rc,2)
    x1=x*r+Rc(1,ii);
    y1=y*r+Rc(2,ii);
    z1=z*r+Rc(3,ii);
    surf(x1,y1,z1,'EdgeColor','none','FaceAlpha',0.2);
    hold on
end
end
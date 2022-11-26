n=100;
p=point_cell{n};
k = convhull(p,'Simplify',true);
show_cylinder(Rc(:,n),Ori(:,n))
trisurf(k,p(:,1),p(:,2),p(:,3),'FaceColor','cyan','FaceAlpha',0.5)

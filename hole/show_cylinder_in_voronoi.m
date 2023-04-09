n=100;
near_pore_cylinder=pore(n).near_cylinder;
for cc=near_pore_cylinder
    p=point_cell{cc};
    k = delaunay(p);
    %     show_cylinder(Rc(:,cc),Ori(:,cc))
    trisurf(k,p(:,1),p(:,2),p(:,3),'FaceColor','cyan','FaceAlpha',0.5)
    hold on
end
show_pore(pore,n,Rc,Ori)

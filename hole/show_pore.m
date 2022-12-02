function show_pore(pore,n,Rc,Ori)
clf
near_pore_cylinder_id=pore(n).near_cylinder;
rc=pore(n).rc;
neighbor=pore(n).real_neighbor;
for ii=1:length(near_pore_cylinder_id)
    show_cylinder(Rc(:,near_pore_cylinder_id(ii)),Ori(:,near_pore_cylinder_id(ii)))
end
plot3(rc(1),rc(2),rc(3),'b*')
cnt=0;
for ii=neighbor
    cnt=cnt+1;
    rc_1=pore(ii).rc;
    plot3(rc_1(1),rc_1(2),rc_1(3),'ko')
    path_tmp=pore(n).path{cnt};
    path_tmp=sortrows(path_tmp')';
    plot3(path_tmp(1,:),path_tmp(2,:),path_tmp(3,:))
end
end
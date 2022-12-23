function show_pore(pore,n,Rc,Ori)
% clf
if pore(n).is_del
    disp('the pore is deleted')
    return
end
near_pore_cylinder_id=pore(n).near_cylinder;
rc=pore(n).rc;
neighbor=pore(n).neighbor;
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
% throat=pore(n).throat;
% plot3(throat(:,1),throat(:,2),throat(:,3),'go')
connective=pore(n).connective;
point=pore(n).area;
for ii=1:length(connective)
    point_tmp=point{ii};
    trisurf(connective{ii},point_tmp(:,1),point_tmp(:,2),point_tmp(:,3),'FaceAlpha',0.5)
end
end
function show_pore(pore_rc,c,v,Rc,Ori,nei_rc,path,throat)
plot3(pore_rc(:,1),pore_rc(:,2),pore_rc(:,3),'*');
hold on
plot3(nei_rc(:,1),nei_rc(:,2),nei_rc(:,3),'go');
show_cylinder(Rc,Ori);
trisurf(c,v(:,1),v(:,2),v(:,3),'EdgeColor','none');
light
for ii=1:length(path{1})
    p=v(path{1}{ii},:);
    plot3(p(:,1),p(:,2),p(:,3))
end
for ii=1:length(throat)
    p=v(throat{ii}(1),:);
    plot3(p(:,1),p(:,2),p(:,3))
end
end
function show_particle(Rc,Ori,v,c,pore_rc,co)
% show_cylinder(Rc,Ori);
hold on
plot3(pore_rc(:,1),pore_rc(:,2),pore_rc(:,3),'*');
map=unique(cell2mat(co));
map(map==0)=[];% 0意味着在内部，不需要上色
color_map=(0:length(map))/length(unique(map));
co=cellfun(@(co1)color_map(1+(1:length(map))*(co1==map')),co,'UniformOutput',false);
co=cellfun(@(c01)c01/max(c01)*255,co,'UniformOutput',false);
for ii=1:length(c)
    trisurf(c{ii},v(:,1),v(:,2),v(:,3),co{ii},'EdgeColor','none');
end
end
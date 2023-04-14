load('D:\毕个业\basic\data\104_plate\basic.mat')
load('D:\毕个业\basic\data\104_plate\pore_data.mat')

[~,idx]=mink(cellfun(@(c)size(c,2),rc_g_cell),100);
%%
% ii=idx(2);
clf
ii=2353;
idx_dcell=idx_dcell_cell{ii};
rc_g=rc_g_cell{ii};
show_cylinder(Rc(:,idx_dcell),Ori(:,idx_dcell),'r')
plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.')
hold on 
pore_near=find(cellfun(@(c)length(intersect(c,idx_dcell))>2,idx_dcell_cell));
for jj=1:length(pore_near)
    if pore_near(jj)==ii
        continue
    end
    rc_g=rc_g_cell{pore_near(jj)};
    idx_dcell=idx_dcell_cell{pore_near(jj)};

    show_cylinder(Rc(:,idx_dcell),Ori(:,idx_dcell),'r')
    plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.','Color',rand(1,3))
end 
load('D:\毕个业\basic\pack\pack104\basic.mat')
load('D:\毕个业\basic\pack\pack104\pore_data.mat')

[~,idx]=mink(cellfun(@(c)size(c,2),rc_g_cell),100);
%%
% ii=idx(2);
clf
ii=5000;
idx_dcell=idx_dcell_cell{ii};
rc_g=rc_g_cell{ii};
[x,y,z]=sphere();
hold on
for ii=1:length(idx_dcell)
    x1=x*r+Rc(1,idx_dcell(ii));
    y1=y*r+Rc(2,idx_dcell(ii));
    z1=z*r+Rc(3,idx_dcell(ii));
    surf(x,y,z,'EdgeColor','none')
    hold on
end
plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.')
hold on
pore_near=find(cellfun(@(c)length(intersect(c,idx_dcell))>2,idx_dcell_cell));
for jj=1:length(pore_near)
    rc_g=rc_g_cell{pore_near(jj)};
    plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.','Color',rand(1,3))
end
axis equal
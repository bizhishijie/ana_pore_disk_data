function h=show_pore(rc_g,r_sphere)
clf
% show_pore(rc_g_cell{1606},r_sphere_cell{1606});
[xs,ys,zs]=sphere(30);
for jj=1:size(rc_g,2)  
    xs2=xs*r_sphere(jj)+rc_g(1,jj);
    ys2=ys*r_sphere(jj)+rc_g(2,jj);
    zs2=zs*r_sphere(jj)+rc_g(3,jj);
    h=surf(xs2,ys2,zs2,ones(size(xs2)));
    set(h,'EdgeAlpha',0)
    hold on
    %         drawnow
end
axis equal
camlight
colormap([0.7 0.7 0.7])
axis off
set(gcf,'color',[1 1 1])
end
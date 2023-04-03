function [] = cylinder_draw_ (rc,ori,D,H,c)
% draw cylinder
% D: diameter
% H: height

%%
[x,y,z]=cylinder(D/2,100);
z=(z-0.5)*H;
rc_grid=[x(:) y(:) z(:)];
V_rot=cross(ori,[0 0 1]');
%dot(ori,[0 0 1]')
th_rot=asin(norm(V_rot)*(2*(dot(ori,[0 0 1]')>0)-1));
M_rot=axang2rotm([V_rot' th_rot]);

rc_grid=rc_grid*M_rot;

rc_grid=rc_grid+repmat(rc',size(rc_grid,1),1);

hold on
K=convhulln(rc_grid);
%h=trisurf(K,rc_grid(:,1),rc_grid(:,2),rc_grid(:,3),ones(size(x))*round(rc(3)));
h=trisurf(K,rc_grid(:,1),rc_grid(:,2),rc_grid(:,3),ones(size(x))*c);
set(h,'EdgeAlpha',0,'FaceAlpha',0.3)

%%
for ii=1:2
    h=plot3(rc_grid(ii:2:end,1),rc_grid(ii:2:end,2),rc_grid(ii:2:end,3),'-');
    set(h,'LineWidth',1.5,'Color',[0.2 0.2 0.2])
end

end


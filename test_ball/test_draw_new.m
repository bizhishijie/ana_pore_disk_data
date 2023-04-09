ii=15;
%% 划分网格
% disp(ii)
idx_tetra=IDX_merge2{ii};
idx_dcell=Td(idx_tetra,:);
idx_dcell=unique(idx_dcell(:));
Rc_tmp=Rc(:,unique(Td_fix(idx_tetra,:)));
%         if any(~is_in(idx_dcell))
%             continue
%         end
Rc_d=ball_sample_point(idx_dcell,:)';
range=[min(Rc_d,[],2) max(Rc_d,[],2)];
[xg,yg,zg]=meshgrid(range(1,1):dx:range(1,2),range(2,1):dx:range(2,2),range(3,1):dx:range(3,2));
rc_g=[xg(:) yg(:) zg(:)]';
% if length(rc_g)>1000000% 超出内存限制，多半是没有合并好
%     continue
% end
clear Rc_d idx_dcell xg yg zg
%% 判断点是否需要留下
idx_out_particle=(all(dist(rc_g',Rc)>r,2))';
idx_in_tetra_loop=false(length(idx_tetra),size(rc_g,2));

parfor jj=1:length(idx_tetra)
    Rc_tetra=ball_sample_point(Td(idx_tetra(jj),:),:)';
    A=Rc_tetra(:,1)-Rc_tetra(:,4);
    B=Rc_tetra(:,2)-Rc_tetra(:,4);
    C=Rc_tetra(:,3)-Rc_tetra(:,4);
    rc_g2=rc_g-Rc_tetra(:,4);

    M3=[A B C]'*[A B C];
    N3=[A B C]'*rc_g2;
    cf3=M3\N3;
    idx_in_tetra_loop(jj,:)=sum(cf3,1)>=0&sum(cf3,1)<=1&all(cf3>0,1);
end
clear A B C M3 N3 cf3
idx_in_tetra=any(idx_in_tetra_loop,1);
idx_keep=idx_out_particle&idx_in_tetra;
rc_g=rc_g(:,idx_keep);
idx_dcell=unique(Td_fix(idx_tetra,:));

clf
plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.')
hold on
for ii=1:length(idx_dcell)
    [x,y,z]=sphere();
    x=x*r+Rc(1,idx_dcell(ii));
    y=y*r+Rc(2,idx_dcell(ii));
    z=z*r+Rc(3,idx_dcell(ii));
    surf(x,y,z,'EdgeColor','none','FaceAlpha',0.1);
end
axis equal
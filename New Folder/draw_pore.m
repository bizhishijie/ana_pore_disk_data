ii=2234;
figure(1);clf;
idx_show=IDX_merge2{ii};
show_cylinder(Rc0(:,unique(Td_fix(idx_show,:))),Ori(:,unique(Td_fix(idx_show,:))),'r')
% for ii=1:length(idx_show)
%     if ~isempty(IDX_merge{idx_show(ii)})
%         p_tmp=[rs_pore(:,IDX_merge{idx_show(ii)}) rs_pore(:,idx_show(ii))];
%         plot3(p_tmp(1,:),p_tmp(2,:),p_tmp(3,:),'b-o')
%         hold on
%     else
%         p_tmp=[rs_pore(:,idx_show(ii))];
%         plot3(p_tmp(1,:),p_tmp(2,:),p_tmp(3,:),'o','MarkerFaceColor','r')
%     end
% end
xlabel('x');ylabel('y');zlabel('z')
axis equal

idx_tetra=IDX_merge2{ii};
idx_dcell=Td(idx_tetra,:);
idx_dcell=unique(idx_dcell(:));
Rc0_tmp=Rc0(:,unique(Td_fix(idx_tetra,:)));
Ori_tmp=Ori(:,unique(Td_fix(idx_tetra,:)));
Rc_d=Rc(:,idx_dcell);
range=[min(Rc_d,[],2) max(Rc_d,[],2)];
[xg,yg,zg]=meshgrid(range(1,1):dx:range(1,2),range(2,1):dx:range(2,2),range(3,1):dx:range(3,2));
rc_g=[xg(:) yg(:) zg(:)]';
clear Rc_d idx_dcell xg yg zg

d1=dot((reshape(rc_g,3,1,[])-Rc0_tmp),repmat(Ori_tmp,1,1,size(rc_g,2)));
d2=vecnorm((reshape(rc_g,3,1,[])-Rc0_tmp)-d1.*repmat(Ori_tmp,1,1,size(rc_g,2)));
d1=max(abs(d1)-h/2,0);d2=max(d2-r,0);
idx_out_particle=permute(all(d1>0|d2>0,[1,2]),[2,3,1]);% 在同一个圆盘的两个面的外面
clear d1 d2
idx_in_tetra_loop=false(length(idx_tetra),size(rc_g,2));
parfor jj=1:length(idx_tetra)
    Rc_tetra=Rc(:,Td(idx_tetra(jj),:));
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

plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.')
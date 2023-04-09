clear
load('basic_418.mat','Rc','idx_eff')

[~,min_idx]=mink(vecnorm(Rc-mean(Rc,2)),100);
Rc=Rc(:,min_idx);
Rc=Rc(:,[
     1
     2
     4
     8]);
D=23.7;
dx=D/20;

[Vv,Cv]=voronoin(Rc');
Td=delaunayn(Rc');

rs_pore=zeros(3,size(Td,1));
dist_pore=zeros(1,size(Td,1));
Nd_idx=cell(1,size(Td,1));
for ii=1:size(Td,1)
    freq=tabulate(cell2mat(Cv(Td(ii,:))'));
    idx=freq(freq(:,2)==4,1);
    idx=idx(idx~=1);% to remove inf
    rs_pore_tmp=Vv(idx,:)';
    rs_pore(:,ii)=rs_pore_tmp;
    dist_pore(ii)=norm(Rc(:,Td(ii,1))-rs_pore_tmp)-D/2;

    Nd_idx{ii}=find(sum(ismember(Td,Td(ii,:)),2)==3);
end

IDX_merge=cell(1,size(Td,1));
for ii=1:size(Td,1)
    n_idx=Nd_idx{ii}';
    dist_p2p=sqrt(sum((rs_pore(:,n_idx)-rs_pore(:,ii)).^2,1));

    IDX_merge{ii}=[n_idx(dist_p2p<dist_pore(ii)&dist_p2p<dist_pore(n_idx)) ii];
end

IDX_merge2 = cell_combine_(IDX_merge);
clear Cv Vv rs_pore dist_pore Nd_idx IDX_merge
%%
[xs,ys,zs]=sphere(80);
% h=surf(xs,ys,zs,ones(size(xs)));
%%
for ii=1:length(IDX_merge2)
% for ii=1
    idx_tetra=IDX_merge2{ii};
    idx_dcell=Td(idx_tetra,:);
    idx_dcell=unique(idx_dcell(:));
    if any(~ismember(idx_dcell,idx_eff))
        continue
    end

    Rc_d=Rc(:,idx_dcell);
    range=[min(Rc_d,[],2) max(Rc_d,[],2)];
    [xg,yg,zg]=meshgrid(range(1,1):dx:range(1,2),range(2,1):dx:range(2,2),range(3,1):dx:range(3,2));
    rc_g=[xg(:) yg(:) zg(:)]';
    
    clear Rc_d range xg yg zg

    idx_out_particle=min(pdist2(rc_g',Rc'),[],2)'>D/2;
    idx_in_tetra_loop=zeros(length(idx_tetra),size(rc_g,2));
    for jj=1:length(idx_tetra)
        Rc_tetra=Rc(:,Td(idx_tetra(jj),:));
        A=Rc_tetra(:,1)-Rc_tetra(:,4);
        B=Rc_tetra(:,2)-Rc_tetra(:,4);
        C=Rc_tetra(:,3)-Rc_tetra(:,4);
        rc_g2=rc_g-Rc_tetra(:,4);

        M3=[A B C]'*[A B C];
        N3=[A B C]'*rc_g2;
        cf3=M3\N3;
        idx_in_tetra_loop(jj,:)=sum(cf3,1)>=0&sum(cf3,1)<=1&all(cf3>0,1)&all(cf3<1,1);
    end
    clear A B C M3 N3 cf3
    idx_in_tetra=any(idx_in_tetra_loop,1);
    idx_keep=idx_out_particle&idx_in_tetra;
    rc_g=rc_g(:,idx_keep);
    plot3(rc_g(1,:),rc_g(2,:),rc_g(3,:),'.');axis equal;hold on
    show_sphere(Rc(:,idx_dcell))
    clear idx_out_particle idx_in_tetra idx_keep

    d_tetra_loop=inf(4*length(idx_tetra),size(rc_g,2));
    idx_surf=nchoosek(unique(reshape(Td(idx_tetra,:),1,[])),3);
    idx_tetra_tmp=Td(idx_tetra,:);
    surf_belong=cell(1,size(idx_surf,1));
    for jj=1:size(idx_surf,1)
        surf_belong{jj}=find(sum(reshape(idx_surf(jj,:),1,1,[])==idx_tetra_tmp,[3,2])==3);
    end
    surf_common=idx_surf(cellfun(@length,surf_belong)==2,:);
   
    for jj=1:length(idx_tetra)
        Rc_tetra=Rc(:,Td(idx_tetra(jj),:));
        surf_common2=surf_common(all(sum(reshape(Td(idx_tetra(jj),:),1,1,[])==surf_common,3),2),:);
        map(Td(idx_tetra(jj),:))=1:4;
        idx_delete=mod(10-sum(map(surf_common2),2)+2,4);% the rows need to be deleted
        idx_delete(idx_delete==0)=4;
        idx_permutation=[1 2 3 4;2 3 4 1;3 4 1 2;4 1 2 3];
        idx_permutation(idx_delete,:)=[];
        for kk=1:size(idx_permutation,1)
            A=Rc_tetra(:,idx_permutation(kk,1))-Rc_tetra(:,idx_permutation(kk,4));
            B=Rc_tetra(:,idx_permutation(kk,2))-Rc_tetra(:,idx_permutation(kk,4));
            C=Rc_tetra(:,idx_permutation(kk,3))-Rc_tetra(:,idx_permutation(kk,4));
            rc_g2=rc_g-Rc_tetra(:,idx_permutation(kk,4));
            M3=[A B C]'*[A B C];
            N3=[A B C]'*rc_g2;
            cf3=M3\N3;
            d_AB=abs(cf3(3,:)*dot(cross(A,B),C)/norm(cross(A,B)));% distacne to ABD
            M2=[A B]'*[A B];
            N2=[A B]'*rc_g2;
            cf2=M2\N2;
            d_AB(sum(cf2,1)<0|sum(cf2,1)>1)=inf;% distance to surface C
            d_tetra_loop((jj-1)*4+kk,:)=d_AB;
        end
    end
   
    d_tetra=min(d_tetra_loop,[],1);
    d_particle=min(pdist2(rc_g',Rc'),[],2)'-D/2;
    r_sphere=min(d_tetra,d_particle);
    clear d_tetra_loop d_AB A B C rc_g2 M3 N3 M2 N2 cf2

    %% remove unnecessary spheres
%     is_rm=false(1,length(rc_g));
%     is_ne=false(1,length(rc_g));
%     while sum(is_rm)/length(is_rm)<0.9
%         disp(sum(is_rm)/length(is_rm));
%         [~,idx_max]=max(r_sphere.*~is_rm);
%         volume_intersect=volume_spherical_cap_(r_sphere(idx_max),r_sphere.*~is_rm,dist(rc_g(:,idx_max)',rc_g.*~is_rm));
%         volume_intersect_rate=volume_intersect./(pi*4/3*r_sphere.^3);
%         is_rm(volume_intersect_rate>0.8)=true;
%         is_rm(idx_max)=true;
%         is_ne(idx_max)=true;
%     end
%     figure(2);clf;hold on
%     r_sphere1=r_sphere(is_ne);
%     rc_g1=rc_g(:,is_ne);
%     for jj=1:size(rc_g1,2)
%         xs2=xs*r_sphere1(jj)+rc_g1(1,jj);
%         ys2=ys*r_sphere1(jj)+rc_g1(2,jj);
%         zs2=zs*r_sphere1(jj)+rc_g1(3,jj);
%         h=surf(xs2,ys2,zs2,ones(size(xs2)));
%         set(h,'EdgeAlpha',0)
% %         drawnow
%     end
%     axis equal
%     camlight

    %%
    [~,idx_max]=max(r_sphere);
    idx_keep=idx_max;
    idx_remain=1:size(rc_g,2);
    while min(r_sphere(idx_keep))/dx>1
        [~,idx_max]=max(r_sphere(idx_remain));
        idx_keep=[idx_keep idx_remain(idx_max)];
        dist_s=sqrt(sum((rc_g(:,idx_remain)-rc_g(:,idx_remain(idx_max))).^2,1));
        idx_remove1=idx_remain(dist_s<=r_sphere(idx_remain(idx_max))*0.7);
        idx_remove2=idx_remain(dist_s<=r_sphere(idx_remain(idx_max))*1.1&...
            r_sphere(idx_remain)<r_sphere(idx_remain(idx_max))*0.3);
        idx_remove=unique([idx_remove1 idx_remove2 idx_remain(idx_max)]);
        idx_remain=idx_remain(~ismember(idx_remain,idx_remove));
%         whos idx_remain idx_remove idx_keep
    end
    LOGIC_draw=false(1,size(rc_g,2));
    LOGIC_draw(idx_keep)=1;

    %%
    figure(2);clf;hold on
    r_sphere1=r_sphere(LOGIC_draw);
    rc_g1=rc_g(:,LOGIC_draw);
    for jj=1:size(rc_g1,2)
        xs2=xs*r_sphere1(jj)+rc_g1(1,jj);
        ys2=ys*r_sphere1(jj)+rc_g1(2,jj);
        zs2=zs*r_sphere1(jj)+rc_g1(3,jj);
        h=surf(xs2,ys2,zs2,ones(size(xs2)));
        set(h,'EdgeAlpha',0)
%         drawnow
    end
    axis equal
    camlight
    colormap([0.7 0.7 0.7])
    axis off
    set(gcf,'color',[1 1 1])
    material([0.5 0.8 0])
    
end




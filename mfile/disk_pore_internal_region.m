%% calculate a proper region for generating random points

clearvars -except pack_num
D=30.75*2/0.8;
H=4.81*2/0.8;
Vp=D^2/4*H*pi;

load_address=['D:\xiachj\research\src disk pore\data\'];
% pack_num=30;

load([load_address num2str(pack_num) '_disk\all_basic_data.mat'],'Rc')

%%
d2=pdist2(Rc',Rc');
d2(d2==0)=inf;
idx_in=find(min(d2,[],1)<D/2);
Rc2=Rc(:,idx_in);
K=convhull(Rc2');

L_cube=max(max(Rc,[],2)-min(Rc,[],2));
Rc_cube=[0 1 1;0 1 -1; 0 -1 1; 0 -1 -1;...
    2 1 1;2 1 -1; 2 -1 1; 2 -1 -1]';
Rc_cube=Rc_cube*L_cube;

Rc_in=Rc2(:,unique(K(:)));

for ii=1:size(K,1)
    Rc_tmp=Rc2(:,K(ii,:));
    n_vec=cross(Rc_tmp(:,2)-Rc_tmp(:,1),Rc_tmp(:,3)-Rc_tmp(:,2));
    n_vec=n_vec/norm(n_vec);
    n_vec=n_vec*sign(dot(n_vec,mean(Rc2,2)-mean(Rc_tmp,2)));
    p_vec=null(n_vec');
    Mrot=[n_vec p_vec];
    Rc_cube2=Mrot*Rc_cube;
    Rc_cube2=Rc_cube2+mean(Rc_tmp,2)+n_vec*D;

    [~,Rc_in]=intersects_polyhedron_(Rc_in,Rc_cube2);
end

%%
d2_in=pdist2(Rc_in',Rc_in');
IDX_near=cell(1,size(Rc_in,2));
for ii=1:size(Rc_in,2)
    IDX_near{ii}=find(d2_in(ii,:)<D/5);
end
IDX_near2=cell_combine_(IDX_near);
Rc_in2=zeros(3,length(IDX_near2));
for ii=1:length(IDX_near2)
    Rc_in2(:,ii)=mean(Rc_in(:,IDX_near2{ii}),2);
end

%%
Rc_in=Rc_in2;
save_address=[load_address num2str(pack_num) '_disk\'];
save([save_address 'internal_polyhedron.mat'],'Rc_in')


%%
figure(1);clf;hold on
K=convhulln(Rc_in');
h=trisurf(K,Rc_in(1,:),Rc_in(2,:),Rc_in(3,:),ones(1,size(Rc_in,2)));
set(h,'FaceAlpha',0.5)
K=convhulln(Rc_in2');
h=trisurf(K,Rc_in2(1,:),Rc_in2(2,:),Rc_in2(3,:),ones(1,size(Rc_in2,2))*2);
set(h,'FaceAlpha',0.5)

drawnow




%% to check
pack_num_loop=31:75;
for ii=pack_num_loop
    load(['D:\xiachj\research\src disk pore\data\' num2str(ii) '_disk\all_basic_data.mat'])
    load(['D:\xiachj\research\src disk pore\data\' num2str(ii) '_disk\internal_polyhedron.mat'])
    idx_eff=[];
    for jj=1:size(Rc,2)
        if inside_polyhedron_(Rc(:,jj),Rc_in)
            idx_eff=[idx_eff jj];
        end
    end
    figure(1);clf;hold on
    plot3(Rc(1,:),Rc(2,:),Rc(3,:),'k.')
    plot3(Rc(1,idx_eff),Rc(2,idx_eff),Rc(3,idx_eff),'ro')
    axis equal
    drawnow
    disp(ii)
    pause
end







%% calculate spherical contact

clearvars -except pack_num file_list
D=30.75*2/0.8;
H=4.81*2/0.8;
Vp=D^2/4*H*pi;

load_address=['D:\毕个业\basic\pack\pack'];
save_address=['D:\毕个业\basic\pack\pack'];
% pack_num=4;

% load([load_address num2str(pack_num) '_disk\all_basic_data.mat'],'Rc','Ori','Vcell')
load([load_address num2str(pack_num) '\basic.mat'],'Rc','Ori')
load([load_address num2str(pack_num) '\edge.mat'],'edge')
Rc_in=edge;

%%
num_point_all=1e5;
disk_pore_contact_spherical_MC_

%%
[p,x]=hist(rs,100);

figure(1);clf
semilogy(x,p/sum(p)/(x(2)-x(1)),'.-')

figure(2);clf
hold on
plot3(rc_p(1,:),rc_p(2,:),rc_p(3,:),'.','MarkerSize',1,'Color',[0.7 0.7 0.7])

rs_thres=16;
idx=find(rs>rs_thres);

plot3(rc_p(1,idx),rc_p(2,idx),rc_p(3,idx),'r.')
grid on
axis equal


%% to remove wrong points due to image processing
d2=pdist2(rc_p(:,idx)',rc_p(:,idx)');
d2(d2==0)=inf;

dis_thres=10;

IDX_near=cell(1,length(idx));
for ii=1:length(idx)
    IDX_near{ii}=find(d2(ii,:)<dis_thres);
end
IDX_pore=cell_combine_ (IDX_near);


pore_size=cellfun('length',IDX_pore);
idx_wrong_point=idx(cell2mat(IDX_pore(pore_size>max(1e2,length(idx)/3))));
idx_wrong_disk=unique(idx_contact(2,idx_wrong_point));%%
idx_correct_point=~ismember(idx_contact(2,:),idx_wrong_disk);
save([save_address num2str(pack_num) '\idx_wrong_disk.mat'],'idx_wrong_disk')
figure(2)
plot3(rc_p(1,idx_wrong_point),rc_p(2,idx_wrong_point),rc_p(3,idx_wrong_point),'k.')

figure(1);hold on
[p,x]=hist(rs(idx_correct_point),100);
semilogy(x,p/sum(p)/(x(2)-x(1)),'.-')


%% calculate precise Phi again
idx_eff=[];
for ii=1:size(Rc,2)
    if inside_polyhedron_(Rc(:,ii),Rc_in)
        idx_eff=[idx_eff ii];
    end
end
idx_eff=idx_eff(~ismember(idx_eff,idx_wrong_disk));
% Phi=length(idx_eff)*Vp/sum(Vcell(idx_eff));


%% large number calculation
num_point_all=1e5;%%
disk_pore_contact_spherical_MC_
idx_correct_point=~ismember(idx_contact(2,:),idx_wrong_disk);

figure(1);hold on
[p,x]=hist(rs(idx_correct_point),100);
semilogy(x,p/sum(p)/(x(2)-x(1)),'.-')

%% save
rc_p=rc_p(:,idx_correct_point);
rs=rs(idx_correct_point);
idx_contact=idx_contact(:,idx_correct_point);

% save_address=[load_address num2str(pack_num) '_disk\'];
% save([save_address 'spherical_contact.mat'],'rs','rc_p','idx_contact','Phi','idx_wrong_disk')




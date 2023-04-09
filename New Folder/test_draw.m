nn=1;

fileList=dir('..\basic\pack\pack*');
load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
load(['..\basic\pack\' fileList(nn).name '\edge.mat'])
load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
load(['..\basic\pack\' fileList(nn).name '\idx_wrong_disk.mat'])
load(['..\basic\pack\' fileList(nn).name '\pore_ball.mat'])

Rc(:,idx_wrong_disk)=[];Ori(:,idx_wrong_disk)=[];is_in(:,idx_wrong_disk)=[];

for ii=1:length(r_sphere_cell)
    clf
    show_pore(rc_g_cell{ii},r_sphere_cell{ii});
    hold on
    particle_id=unique(reshape(pore_particle_id((IDX_merge2{ii}'-1)*4+(1:4)),1,[]));
    show_cylinder(Rc(:,particle_id),Ori(:,particle_id))
end
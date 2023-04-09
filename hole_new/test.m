load('D:\毕个业\basic\pack\pack12\voro.mat')
load('D:\毕个业\basic\pack\pack12\pore_1.mat')
load('D:\毕个业\basic\pack\pack12\volume.mat')
%%
volume(volume<0)=0;
pore_is_in=cellfun(@(c)all(is_in(c)),pore_c_id);
volume_new=volume(pore_is_in'&~is_del);
hist(volume_new)
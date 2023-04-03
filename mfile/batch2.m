file_list=dir('D:\毕个业\basic\pack\pack*');
for ii=1:length(file_list)
    pack_num=file_list(ii).name;
    disp(pack_num);
    pack_num(1:4)=[];
    pack_num=str2num(pack_num);
    disk_pore_contact_spherical
end
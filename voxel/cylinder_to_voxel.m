function voxel=cylinder_to_voxel(Rc,Ori,box,cylinder_size,voxel_num)
voxel_num=int32(voxel_num);
box_min=box(:,1);
box_max=box(:,2);
voxel=false(voxel_num,voxel_num,voxel_num);
for ii=int32(1:voxel_num^3)

    vector_3=ii/voxel_num/voxel_num;
    vector_2=(ii-vector_3*voxel_num*voxel_num)/voxel_num;
    vector_1=(ii-vector_3*voxel_num*voxel_num-vector_2*voxel_num);
    vector=double([vector_1;vector_2;vector_3])/double(voxel_num);
    % 以上三行利用了 int 类型的整除
    point_voxel=vector.*(box_max-box_min)+box_min;
    if is_point_in_cylinder(point_voxel,Rc,Ori,cylinder_size)
        voxel(ii)=true;
    end
end
end
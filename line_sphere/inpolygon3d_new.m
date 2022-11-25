function is_in=inpolygon3d_new(point,rc_list)
idx_face=convhulln(rc_list');
rc2_c=mean(rc_list,2);
face_side=zeros(1,size(idx_face,1));
for jj=1:size(idx_face,1)
    n_vec=cross(rc_list(:,idx_face(jj,1))-rc_list(:,idx_face(jj,2)),rc_list(:,idx_face(jj,2))-rc_list(:,idx_face(jj,3)));
    n_vec=n_vec/norm(n_vec);
    face_side_tmp=dot(point-rc_list(:,idx_face(jj,1)),n_vec)*dot(rc2_c-rc_list(:,idx_face(jj,1)),n_vec);
    if face_side_tmp>0
        face_side(jj)=1;
    else
        break
    end
end
is_in=(sum(face_side==0)==0);
% 是1表明在内部
end
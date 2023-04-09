function poly_inside=inside_polyhedron_(rc,rc_p)
% to determine whether a point is inside a polyhedron

K=convhulln(rc_p');
rc_pm=mean(rc_p,2);
face_side=zeros(1,size(K,1));
for jj=1:size(K,1)
    n_vec=cross(rc_p(:,K(jj,1))-rc_p(:,K(jj,2)),rc_p(:,K(jj,2))-rc_p(:,K(jj,3)));
    n_vec=n_vec/norm(n_vec);
    face_side_tmp=dot(rc-rc_p(:,K(jj,1)),n_vec)*dot(rc_pm-rc_p(:,K(jj,1)),n_vec);
    if face_side_tmp>0
        face_side(jj)=1;
    else
        break
    end
end
poly_inside=(sum(face_side==0)==0);% 1 for inside
end
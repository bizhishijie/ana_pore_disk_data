function [V_inter,rc_inter_all] = intersects_polyhedron_(rc1,rc2)
%% intersection between two polyhedron

% rc1=rand(3,10);
% rc2=rand(3,10)+0.2;
%%
K1=convhulln(rc1');
K2=convhulln(rc2');

%%
% figure(1);clf;hold on
% h=trisurf(K1,rc1(1,:),rc1(2,:),rc1(3,:),ones(1,size(rc1,2)));
% set(h,'FaceAlpha',0.1)
% h=trisurf(K2,rc2(1,:),rc2(2,:),rc2(3,:),ones(1,size(rc2,2))*2);
% set(h,'FaceAlpha',0.1)
% axis equal

%%

idx_line1=[K1(:,[1 2]);K1(:,[2 3])];
idx_face1=K1;
idx_line2=[K2(:,[1 2]);K2(:,[2 3])];
idx_face2=K2;

rc_inter_line=[];

for ii=1:size(idx_line1,1)
    for jj=1:size(idx_face2,1)
        rc_inter_tmp=intersects_line_face_ (rc1(:,idx_line1(ii,:)),rc2(:,idx_face2(jj,:)));
        rc_inter_line=[rc_inter_line rc_inter_tmp];
    end
end
for ii=1:size(idx_line2,1)
    for jj=1:size(idx_face1,1)
        rc_inter_tmp=intersects_line_face_ (rc2(:,idx_line2(ii,:)),rc1(:,idx_face1(jj,:)));
        rc_inter_line=[rc_inter_line rc_inter_tmp];
    end
end

%%
rc_inter_corner=[];
rc2_c=mean(rc2,2);
for ii=1:size(rc1,2)
    face_side=zeros(1,size(idx_face2,1));
    for jj=1:size(idx_face2,1)
        n_vec=cross(rc2(:,idx_face2(jj,1))-rc2(:,idx_face2(jj,2)),rc2(:,idx_face2(jj,2))-rc2(:,idx_face2(jj,3)));
        n_vec=n_vec/norm(n_vec);
        face_side_tmp=dot(rc1(:,ii)-rc2(:,idx_face2(jj,1)),n_vec)*dot(rc2_c-rc2(:,idx_face2(jj,1)),n_vec);
        if face_side_tmp>0
            face_side(jj)=1;
        else
            break
        end
    end
    if sum(face_side==0)==0
        rc_inter_corner=[rc_inter_corner rc1(:,ii)];
    end
end
%% 判断一个点是否在内侧
rc1_c=mean(rc1,2);
for ii=1:size(rc2,2)
    face_side=zeros(1,size(idx_face1,1));
    for jj=1:size(idx_face1,1)
        n_vec=cross(rc1(:,idx_face1(jj,1))-rc1(:,idx_face1(jj,2)),rc1(:,idx_face1(jj,2))-rc1(:,idx_face1(jj,3)));
        n_vec=n_vec/norm(n_vec);
        face_side_tmp=dot(rc2(:,ii)-rc1(:,idx_face1(jj,1)),n_vec)*dot(rc1_c-rc1(:,idx_face1(jj,1)),n_vec);
        if face_side_tmp>0
            face_side(jj)=1;
        else
            break
        end
    end
    if sum(face_side==0)==0
        rc_inter_corner=[rc_inter_corner rc2(:,ii)];
    end
end

%%
% plot3(rc_inter_line(1,:),rc_inter_line(2,:),rc_inter_line(3,:),'ro')
% plot3(rc_inter_corner(1,:),rc_inter_corner(2,:),rc_inter_corner(3,:),'mo')

rc_inter_all=[rc_inter_line rc_inter_corner];
if ~isempty(rc_inter_all)
    [~,V_inter]=convhulln(rc_inter_all');
else
    V_inter=0;
end

end
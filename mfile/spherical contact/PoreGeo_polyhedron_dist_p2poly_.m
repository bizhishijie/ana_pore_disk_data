function t = PoreGeo_polyhedron_dist_p2poly_ (rc_p,rc,ori,L,type)
%求空间中一点与多面体之间的最短距离
%rc_p：空间中一点。rc：多面体质心。ori：多面体取向。L：多面体边长。
%%
switch type
    case'cube'
        rc_nod0=[1 1 1;-1 1 1;1 -1 1;1 1 -1;1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]';
        face_point=[1 2 7 3;1 2 6 4;1 3 5 4;2 6 8 7;3 5 8 7;4 5 8 6]';
        idx_edge=[1,2;1,3;1,4;2,6;2,7;3,5;3,7;4,5;4,6;5,8;6,8;7,8;];
        R_p=sqrt(3)/2*L;
        dist_n2c=sqrt(2)/2*L;
        th_sym=pi/2;
    case'octa'
        rc_nod0=[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]';
        face_point=[3,5,2;5,3,1;6,3,2;3,6,1;4,5,1;5,4,2;4,6,2;6,4,1]';
        idx_edge=[1,3;1,4;1,5;1,6;2,3;2,4;2,5;2,6;3,5;3,6;4,5;4,6;];
        R_p=L/sqrt(2);
        dist_n2c=L/sqrt(3);
        th_sym=pi*2/3;
    case'dodeca'
        m=(1+sqrt(5))/2;
        rc_nod0=[0 m 1/m;0 -m 1/m;0 m -1/m;0 -m -1/m;...
            1/m 0 m;-1/m 0 m;1/m 0 -m;-1/m 0 -m;...
            m 1/m 0;-m 1/m 0;m -1/m 0;-m -1/m 0;...
            1 1 1;-1 1 1;1 -1 1;1 1 -1;1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]';
        face_point=[1 13 5 6 14;1 14 10 18 3;1 3 16 9 13;...
            19 2 4 20 12;11 17 4 2 15;5 15 2 19 6;...
            3 18 8 7 16;8 20 4 17 7;13 9 11 15 5;...
            6 19 12 10 14;16 7 17 11 9;10 12 20 8 18]';
        idx_edge=[1,3;1,13;1,14;2,4;2,15;2,19;...
            3,16;3,18;4,17;4,20;5,6;5,13;...
            5,15;6,14;6,19;7,8;7,16;7,17;...
            8,18;8,20;9,11;9,13;9,16;10,12;...
            10,14;10,18;11,15;11,17;12,19;12,20];
        R_p=L*sqrt(3)*(1+sqrt(5))/4;
        dist_n2c=L/2/sin(pi/5);
        th_sym=pi*2/5;
    case'icosa'
        m=sqrt(50-10*sqrt(5))/10;
        n=sqrt(50+10*sqrt(5))/10;
        rc_nod0=[m 0 n;-m 0 n;m 0 -n;-m 0 -n;...
            0 n m;0 -n m;0 n -m;0 -n -m;...
            n m 0;n -m 0;-n m 0;-n -m 0]';
        face_point=[7,3,9;4,7,11;7,4,3;2,1,6;...
            3,10,9;10,1,9;1,10,6;1,5,9;...
            5,7,9;7,5,11;5,2,11;2,5,1;...
            2,12,11;12,4,11;12,2,6;8,10,3;...
            4,8,3;10,8,6;8,12,6;12,8,4]';
        idx_edge=[1,2;1,5;1,6;1,9;1,10;...
            2,5;2,6;2,11;2,12;3,4;...
            3,7;3,8;3,9;3,10;4,7;...
            4,8;4,11;4,12;5,7;5,9;
            5,11;6,8;6,10;6,12;7,9;...
            7,11;8,10;8,12;9,10;11,12];
        R_p=L*sqrt(10+2*sqrt(5))/4;
        dist_n2c=L/sqrt(3);
        th_sym=pi*2/3;
end

%%
Mrot=eul2rotm(ori','ZYZ');
rc_nod=rc_nod0/norm(rc_nod0(:,1))*R_p;
rc_nod=Mrot*rc_nod+repmat(rc,1,size(rc_nod,2));

%% plane
nor_vec=cross(rc_nod(:,face_point(2,:))-rc_nod(:,face_point(1,:)),...
    rc_nod(:,face_point(3,:))-rc_nod(:,face_point(2,:)));
nor_vec=nor_vec./repmat(sqrt(sum(nor_vec.^2,1)),3,1);
rc_corner=rc_nod(:,face_point(1,:));
rc_face_center=zeros(3,size(face_point,2));
for ii=1:size(face_point,2)
    rc_face_center(:,ii)=mean(rc_nod(:,face_point(:,ii)),2);
end

t_tmp1=dot(rc_corner-rc_p,nor_vec);

rc_inter=rc_p+nor_vec.*t_tmp1;

vec_t2c=rc_inter-rc_face_center;
dist_t2c=sqrt(sum(vec_t2c.^2,1));
vec_t2c=vec_t2c./repmat(dist_t2c,3,1);
vec_n2c=rc_corner-rc_face_center;
vec_n2c=vec_n2c/dist_n2c;

th=acos(dot(vec_t2c,vec_n2c));
th_true=~isreal(th);
if sum(th_true)>0
    t1=[1;1];
else
    th=rem(th+th_sym/2,th_sym);
    th=min(th,th_sym-th);

    idx_inter1=dist_t2c<=dist_n2c*cos(th_sym/2)./cos(th);
    t1=t_tmp1(idx_inter1);

end
[t1,~]=sort(abs(t1));

if ~isempty(t1)
    t1=t1(1);

else
    t1=0;
end


%% line
idx_edge=idx_edge';
edge_nod=reshape(rc_nod(:,idx_edge(:)),6,[]);


t_tmp2=dot((rc_p-edge_nod(1:3,:)),(edge_nod(4:6,:)-edge_nod(1:3,:)))./L;
ver_point=edge_nod(1:3,:)+t_tmp2.*(edge_nod(4:6,:)-edge_nod(1:3,:));
t2=(sum((rc_p-ver_point).^2,1)).^(1/2);

idx_inter2=(t_tmp2>0) & (t_tmp2<L);

[t2,~]=sort(t2(idx_inter2));
if ~isempty(t2)
    t2=t2(1);
else
    t2=0;
end

%% point
[t3,~]=sort((sum((rc_nod-rc_p).^2,1)).^(1/2));
t3=t3(1);


%%
t_all=[t1,t2,t3];

t_all=t_all(t_all>0);
t=sort(t_all);
t=t(1);

end
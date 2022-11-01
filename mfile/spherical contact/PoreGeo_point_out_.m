function t = PoreGeo_point_out_ (rc_p,rc,ori,L,type)
%% 判断空间中一点是否在多面体外部，0表示在内部，1表示在外部。
%%
switch type
    case 'cube'
        rc_nod0=[1 1 1;-1 1 1;1 -1 1;1 1 -1;1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]';
        face_point=[1 2 7 3;1 2 6 4;1 3 5 4;2 6 8 7;3 5 8 7;4 5 8 6]';
        R_p=sqrt(3)/2*L;
        dist_n2c=sqrt(2)/2*L;
        th_sym=pi/2;
    case'octa'
        rc_nod0=[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]';
        face_point=[3,5,2;5,3,1;6,3,2;3,6,1;4,5,1;5,4,2;4,6,2;6,4,1]';
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
        R_p=L*sqrt(10+2*sqrt(5))/4;
        dist_n2c=L/sqrt(3);
        th_sym=pi*2/3;
end

%%
n_vec=[0;0;1];
Mrot=eul2rotm(ori','ZYZ');
rc_nod=rc_nod0/norm(rc_nod0(:,1))*R_p;
rc_nod=Mrot*rc_nod+repmat(rc,1,size(rc_nod,2));

nor_vec=cross(rc_nod(:,face_point(2,:))-rc_nod(:,face_point(1,:)),...
    rc_nod(:,face_point(3,:))-rc_nod(:,face_point(2,:)));
nor_vec=nor_vec./repmat(sqrt(sum(nor_vec.^2,1)),3,1);
rc_corner=rc_nod(:,face_point(1,:));
rc_face_center=zeros(3,size(face_point,2));
for ii=1:size(face_point,2)
    rc_face_center(:,ii)=mean(rc_nod(:,face_point(:,ii)),2);
end

t_tmp=dot(rc_corner-rc_p,nor_vec)./...
    dot(repmat(n_vec,1,size(face_point,2)),nor_vec);

rc_inter=rc_p+n_vec*t_tmp;

vec_t2c=rc_inter-rc_face_center;
dist_t2c=sqrt(sum(vec_t2c.^2,1));
vec_t2c=vec_t2c./repmat(dist_t2c,3,1);
vec_n2c=rc_corner-rc_face_center;
vec_n2c=vec_n2c/dist_n2c;

th=acos(dot(vec_t2c,vec_n2c));
th_true=~isreal(th);
if sum(th_true)>0
    t=0;
else
    th=rem(th+th_sym/2,th_sym);
    th=min(th,th_sym-th);

    idx_inter=dist_t2c<=dist_n2c*cos(th_sym/2)./cos(th);

    if sum(idx_inter)==2
        t=t_tmp(idx_inter);
        t=t(1)*t(2)>0;
        
%         plot3(rc_p(1),rc_p(2),rc_p(3),'o')
%         point1=rc_p+t_tmp(idx_inter).*n_vec;
%         plot3(point1(1,:),point1(2,:),point1(3,:),'*')
%         point2=rc_p+t_tmp(idx_inter).*n_vec*3;
%         plot3(point2(1,:),point2(2,:),point2(3,:),'LineWidth',2)
    else
        t=1;
    end
end
end


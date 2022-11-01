function [dist_s]=PoreGeo_rs_poly_(Rc,Ori,L,resolution,type,N_test)
%% 计算正多面体packing的rs
% dist_s孔隙中测试点到颗粒的最短距离，单位为等体积球直径
% L多面体边长，N_test测试点数目，resolution分辨率（mm/pxiel）
%%
switch type
    case'cube'
        ri=1/2;%内切球半径(颗粒边长为一)
        ru=sqrt(3)/2;%外切球半径
        r_ev=(3/2)^(1/2);%等体积球半径

    case'octa'
        ri=6^(1/2)/6;
        ru=2^(1/2)/2;
        r_ev=(2^(1/2)/4/pi)^(1/3);
    case'dodeca'
        ri=1/2*(5/2+11/10*5^(1/2))^(1/2);
        ru=3^(1/2)/4*(1+5^(1/2));
        r_ev=nthroot(3/16/pi*(15+7*sqrt(5)),3);

    case'icosa'
        ri=(3*3^(1/2)+15^(1/2))/12;
        ru=(10+2*5^(1/2))^(1/2)/4;
        r_ev=((15+5*5^(1/2))/16/pi)^(1/3);
end


%% modify coordinates（xz的数据需要调整方向）
rot_dir=[0 1 0];
rot_angle=1.8/180*pi;
s=cos(rot_angle/2);
t=sin(rot_angle/2);
rot_dir=rot_dir/norm(rot_dir);
x=rot_dir(1)*t;y=rot_dir(2)*t;z=rot_dir(3)*t;
rot_M=[s^2+x^2-y^2-z^2 2*x*y-2*s*z 2*x*z+2*s*y;...
    2*x*y+2*s*z s^2-x^2+y^2-z^2 2*y*z-2*s*x;...
    2*x*z-2*s*y 2*y*z+2*s*x s^2+z^2-x^2-y^2];

% resolution=1;%mm/pxiel
Rc=Rc/(L/resolution);% set edge length = 1

Rc=Rc-mean(Rc,2);
Rc=rot_M*Rc;


%%
Rc_range=[min(Rc,[],2),max(Rc,[],2)];
dis_bdr=min([Rc-repmat(Rc_range(:,1),1,size(Rc,2));...
    repmat(Rc_range(:,2),1,size(Rc,2))-Rc],[],1);

dis_bdr=dis_bdr/r_ev;
idx_eff=find(dis_bdr>=3);

Rc_range2=[min(Rc(:,idx_eff),[],2),max(Rc(:,idx_eff),[],2)];


%%

dist_s=zeros(1,N_test);
Rc_test=zeros(3,N_test);
N_pore=0;

tic
while N_pore<N_test
    rc_test=rand(3,1).*(Rc_range2(:,2)-Rc_range2(:,1))+Rc_range2(:,1);
    dist_p2poly=sqrt(sum((Rc-rc_test).^2,1));
    idx_close1=find(dist_p2poly>ri & dist_p2poly<ru);


    if min(dist_p2poly)>ri
        t=1;  
        for ii=1:length(idx_close1)
            t=t && PoreGeo_point_out_(rc_test,Rc(:,idx_close1(ii)),Ori(:,idx_close1(ii)),1,type);

        end
        
        if t
            idx_close2=find(dist_p2poly<(2*r_ev));

            if isempty(idx_close2)
                [~,idx_close2]=sort(dist_p2poly);
                idx_close2=idx_close2(1);
            end
            dist_tmp=ones(1,length(idx_close2));
                for jj=1:length(idx_close2)
                    dist_tmp(jj)=PoreGeo_polyhedron_dist_p2poly_ (rc_test,Rc(:,idx_close2(jj)),Ori(:,idx_close2(jj)),1,type);
                end
            dist_tmp=min(dist_tmp);
            N_pore=N_pore+1;
            dist_s(N_pore)=dist_tmp/(2*r_ev);
            Rc_test(:,N_pore)=rc_test;
        end

    end
end
toc
end
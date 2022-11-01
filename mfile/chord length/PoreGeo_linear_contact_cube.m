clear;clc
tic
load_address=('.\cube_1_1\pack23\');
file_name=('basic.mat');

load([load_address file_name])


%% modify coordinates
rot_dir=[0 1 0];
rot_angle=1.8/180*pi;
s=cos(rot_angle/2);
t=sin(rot_angle/2);
rot_dir=rot_dir/norm(rot_dir);
x=rot_dir(1)*t;y=rot_dir(2)*t;z=rot_dir(3)*t;
rot_M=[s^2+x^2-y^2-z^2 2*x*y-2*s*z 2*x*z+2*s*y;...
    2*x*y+2*s*z s^2-x^2+y^2-z^2 2*y*z-2*s*x;...
    2*x*z-2*s*y 2*y*z+2*s*x s^2+z^2-x^2-y^2];

L=30;% mm
resolution=1.1;%mm/pxiel
Rc=Rc/(L/resolution);% set cube edge length = 1

Rc=Rc-mean(Rc,2);
Rc=rot_M*Rc;


%% boundary
vessel_r=mean(max(Rc(1:2,:),[],2)-min(Rc(1:2,:),[],2))/2;
vessel_h=[min(Rc(3,:));max(Rc(3,:))];

dis_bdr=[vessel_r-sqrt(sum(Rc(1:2,:).^2,1));...
    Rc(3,:)-vessel_h(1,:);vessel_h(2,:)-Rc(3,:)];

idx_eff=find(min(dis_bdr,[],1)>=0.75);
Rc_eff=Rc(:,idx_eff);
Ori_eff=Ori(:,idx_eff);
Rc_range2=[min(Rc(:,idx_eff),[],2),max(Rc(:,idx_eff),[],2)];


%%
N_test=1e3;
num_th=4;
dist_c=cell(num_th,0);
dist_l=zeros(num_th,N_test);
R_equal_volume=nthroot(3/4/pi,3);

for ii=1:num_th
    N_pore=0;
    N_chord=1;
    while N_pore<N_test
        rc_test=rand(3,1).*(Rc_range2(:,2)-Rc_range2(:,1))+Rc_range2(:,1);
        ph=2*pi*rand(1);
        th=acos(rand(1)/num_th+(ii-1)/num_th);
        n_vec=[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)]';

        if norm(rc_test(1:2))<=vessel_r-1
            dist_lateral=sqrt(sum(cross(rc_test-Rc,...
                repmat(n_vec,1,size(Rc,2))).^2,1));
            idx_eff_close=find(dist_lateral<=1);

            t_loop=zeros(2,length(idx_eff_close));
            for jj=1:length(idx_eff_close)
                t_loop(:,jj)=PoreGeo_polyhedron_intersect_(rc_test,n_vec,...
                    Rc(:,idx_eff_close(jj)),Ori(:,idx_eff_close(jj)),1,'cube');
            end

            t_loop=t_loop(:,t_loop(1,:)~=0);
            t_loop=sort(t_loop,1);
            [~,idx_t]=sort(t_loop(1,:),2);
            t_loop=t_loop(:,idx_t);

            %排除颗粒重叠带来的影响
            t_loop2=t_loop;
            for jj=2:size(t_loop,2)
                t2_max=max(t_loop(2,1:jj-1));
                if t_loop(1,jj)<t2_max
                    t_loop(1,jj)=-inf;
                    t_loop(2,jj)=max(t2_max,t_loop(2,jj));
                    t_loop(2,jj-1)=-inf;
                end
            end
            t_loop=t_loop(:);
            t_loop=reshape(t_loop(~isinf(t_loop)==1),2,[]);


            dist_c{ii,N_chord}=(t_loop(1,2:end)-t_loop(2,1:end-1))./(R_equal_volume*2);
            N_chord=N_chord+1;

            if sum(t_loop(1,:).*t_loop(2,:)<0)==0&&~isempty(find(t_loop>0,1))
                N_pore=N_pore+1;
                dist_l(ii,N_pore)=(min(t_loop(t_loop>0)))./(R_equal_volume*2);
            end
        end
    end
end
toc



%%
figure(1);clf
hold on
th_num=char('cos(th) [0,0.25] 水平方向','cos(th) [0.25,0.50]',...
    'cos(th) [0.50,0.75]','cos(th) [0.75,1] 竖直方向');
for ii=1:4
    dist_chord=cell2mat(dist_c(ii,:));
    [pc,xc]=hist(dist_chord,100);
    pc=pc/sum(pc)/(xc(2)-xc(1));
    plot(xc,pc,'Linewidth',2,'DisplayName',th_num(ii,:))
end

% for ii=1:5
%     [pl,xl]=hist(dist_l(ii,:),100);
%     pl=pl/sum(pl)/(xl(2)-xl(1));
%     plot(xl,pl,color_idx(ii),'Linewidth',2)
% end

set(gca,'yscale','log')
xlabel('chord length/D')
ylabel('P(chord)')
title('Cube    P(chord_length)     Phi=0.750')
legend

%%

% save_address=('C:\Users\WCY\Desktop\绘图数据\ch\');
% file_name=('chord_cube_2.mat');
%
% save([save_address file_name],'dist_c')







%% 3D image
% figure(2);clf
% hold on
% 
% for ii=1:length(idx_eff_close)
%     cube_draw_(Rc_eff(:,idx_eff_close(ii)),Ori_eff(:,idx_eff_close(ii)),1,2,0.5)
% end
% 
% rc_line=[rc_test-10*n_vec,rc_test+10*n_vec]';
% sl=line(rc_line(:,1),rc_line(:,2),rc_line(:,3));
% sl.LineWidth=2;
% sl.Color='red';
% 
% t_sl=t_loop(t_loop~=0);
% rc_sl=rc_test+t_sl'.*n_vec;
% plot3(rc_sl(1,:),rc_sl(2,:),rc_sl(3,:),'o');
% 
% view(3)
% axis equal


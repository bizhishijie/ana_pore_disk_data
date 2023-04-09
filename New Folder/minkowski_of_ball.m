% clear
% load('D:\毕个业\basic\pack\pack104\pore_ball.mat')
% load('pore_ball.mat')
% show_pore(rc_g_cell{2},r_sphere_cell{2})
v_pore_list=zeros(1,length(r_sphere_cell));
eccentricity_list=zeros(1,length(r_sphere_cell));
%% 蒙特卡洛法
% for ii=1:length(r_sphere_cell)
%     disp(ii)
%     r_sphere=r_sphere_cell{ii};
%     range=range_cell{ii};
%     rc_g=rc_g_cell{ii};
%     rand_num=1e4;
%     rand_p=rand(3,rand_num).*(range(:,2)-range(:,1))+range(:,1);
%     id_in=any(dist(rand_p',rc_g)<r_sphere,2);
%     v_pore=sum(id_in)/rand_num*prod(range(:,2)-range(:,1));
%     p_in=rand_p(:,id_in);% 每个点的相对坐标
%     p_in=p_in-mean(p_in,2);
%     mk=diag(p_in*p_in');
%     eccentricity_list(ii)=abs(min(mk)/max(mk));
%     v_pore_list(ii)=v_pore;
% end
%%
h=4.81*2/0.8;
dx=h/20;
for ii=1:length(rc_g_cell)
    rc_g_tmp=rc_g_cell{ii};
    v_pore=dx^3*length(rc_g_tmp);

    mk=diag(rc_g_tmp*rc_g_tmp');
    v_pore_list(ii)=v_pore;
    eccentricity_list(ii)=abs(min(mk)/max(mk));
end
figure(1)
hist_num=100;
x=linspace(0,6e4,hist_num);
y=hist(v_pore_list,x);
plot(x,y/sum(y))
set(gca,'YScale','log')
% loglog(x,y/sum(y))
set(gca,'YScale','log')
figure(2)
hist_num=20;
x=linspace(0,1,hist_num);
y=hist(eccentricity_list,hist_num);
plot(x,y)
set(gca,'YScale','log')
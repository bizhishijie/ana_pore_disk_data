clear;clc
load_address=['E:\research\analysis\sphere\sphere_date_xz\'];
file_name=['sphere_37'];
load([load_address file_name '.mat'])


%%
D=15;
Rc_range=[min(Rc_p,[],2),max(Rc_p,[],2)];
dis_bdr=min([Rc_p-repmat(Rc_range(:,1),1,size(Rc_p,2));...
    repmat(Rc_range(:,2),1,size(Rc_p,2))-Rc_p],[],1);
dis_bdr=dis_bdr/D;
idx_eff=find(dis_bdr>=1.5);
Rc_p_eff=Rc_p(:,idx_eff);
x
Rc_range2=[min(Rc_p(:,idx_eff),[],2),max(Rc_p(:,idx_eff),[],2)];


%%
N_test=1e5;
num_th=4;
dist_c=cell(num_th,0);
dist_l=zeros(num_th,N_test);

for ii=1:num_th
    N_pore=0;
    N_chord=1;   
while N_pore<N_test
    rc_test=rand(3,1).*(Rc_range2(:,2)-Rc_range2(:,1))+Rc_range2(:,1);
    ph=2*pi*rand(1);
    th=acos(rand(1)/num_th+(ii-1)/num_th);
    n_vec=[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)]';
    
    dist_lateral=sqrt(sum(cross(Rc_p-rc_test,repmat(n_vec,1,size(Rc_p,2))).^2,1));
    
    idx_close=find(dist_lateral<=D/2);
    t_loop=zeros(1,2*length(idx_close));
    for jj=1:length(idx_close)
        B=2*dot(Rc_p(:,idx_close(jj))-rc_test,n_vec);
        C=norm(Rc_p(:,idx_close(jj))-rc_test)^2-(D/2)^2;
        t=[(-B+sqrt(B^2-4*C))/2 (-B-sqrt(B^2-4*C))/2];
        t_loop(2*jj-1:2*jj)=t;
    end
    
    t_loop=flip(t_loop);%%
    t_loop=reshape(t_loop,2,length(idx_close));%%
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
    
    
%     t_loop=reshape(t_loop,2,length(idx_close));
%     t_loop=sort(t_loop,1);
%     [~,idx_s]=sort(t_loop(1,:));
%     t_loop=t_loop(:,idx_s);
    
    chord=t_loop(1,2:end)-t_loop(2,1:end-1);
    dist_c{ii,N_chord}=chord./D;
    
    N_chord=N_chord+1;
    if sum(t_loop(1,:).*t_loop(2,:)<0)==0&&~isempty(find(t_loop>0,1))
        N_pore=N_pore+1;
        dist_l(ii,N_pore)=min(t_loop(t_loop>0))./D;
    end
end
end


figure(1);clf
hold on

th_num=char('cos(th) [0,0.25]','cos(th) [0.25,0.50]','cos(th) [0.50,0.75]','cos(th) [0.75,1]');
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
title('Sphere        P（chord length）')
legend


%%
% save_address=('C:\Users\WCY\Desktop\绘图数据\ch\');
% file_name=('chord_sphere.mat');
% 
% save([save_address file_name],'dist_c')


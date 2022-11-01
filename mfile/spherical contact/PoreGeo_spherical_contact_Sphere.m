clear;clc
load_address=['.\'];
file_name=['sphere_37'];
load([load_address file_name '.mat'])

%%


D=15;
% D2=ceil(D/2+1);
% [xm,ym,zm]=meshgrid(-D2:D2,-D2:D2,-D2:D2);


Rc_range=[min(Rc_p,[],2),max(Rc_p,[],2)];
dis_bdr=min([Rc_p-repmat(Rc_range(:,1),1,size(Rc_p,2));...
    repmat(Rc_range(:,2),1,size(Rc_p,2))-Rc_p],[],1);
dis_bdr=dis_bdr/D;
idx_eff=find(dis_bdr>=1.5);

Rc_range2=[min(Rc_p(:,idx_eff),[],2),max(Rc_p(:,idx_eff),[],2)];


%%
N_test=1e4;
dist_s=zeros(1,N_test);
N_pore=0;
while N_pore<N_test
    rc_test=rand(3,1).*(Rc_range2(:,2)-Rc_range2(:,1))+Rc_range2(:,1);
    dist_tmp=min(sqrt(sum((Rc_p-rc_test).^2,1)))-D/2;
    if dist_tmp>=0
        N_pore=N_pore+1;
        dist_s(N_pore)=dist_tmp/D;
    end
end


%%
[ps,xs]=hist(dist_s,100);
ps=ps/sum(ps)/(xs(2)-xs(1));
figure
plot(xs,ps,'o-')
set(gca,'yscale','log')
xlabel('rs/D')
ylabel('p(rs)')
title('Sphere Hs(z)')


%%
% save_address=['C:\Users\WCY\Desktop\绘图数据\sc\'];
% file_name=['spherical_sphere'];
% save([save_address file_name '.mat'],'dist_s');

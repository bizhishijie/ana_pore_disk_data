clear
pack_num_loop=[1:29];

Phi_loop=zeros(1,length(pack_num_loop));
rsm_loop=zeros(1,length(pack_num_loop));

for ii=1:length(pack_num_loop)
    load_address=['D:\xiachj\research\src disk pore\data\' num2str(pack_num_loop(ii)) '_disk\'];
    load([load_address 'spherical_contact.mat'],'Phi','rs')
    Phi_loop(ii)=Phi;
    rsm_loop(ii)=mean(rs);
end


%% compaction
clear
pack_num_loop=[1:30];
% pack_num_loop=[48:50 54:57 59 60 65 66 68 70 73:75];% with crystal order

D=30.75*2/0.8;
H=4.81*2/0.8;
Vp=D^2/4*H*pi;

Phi_loop=zeros(8,length(pack_num_loop));
rsm_loop=zeros(8,length(pack_num_loop));
s2_loop=zeros(8,length(pack_num_loop));

Phi_range=linspace(0.525,0.73,10);
% Phi_range=linspace(0.59,0.75,7);
rs_ALL=cell(1,length(Phi_range)-1);
vcell_ALL=cell(1,length(Phi_range)-1);

for ii=1:length(pack_num_loop)
    load_address=['D:\xiachj\research\src disk pore\data\' num2str(pack_num_loop(ii)) '_disk\'];
    load([load_address 'spherical_contact.mat'],'Phi','rs','idx_wrong_disk','rc_p','idx_contact')
    load([load_address 'all_basic_data.mat'],'Rc','Ori','Vcell')
    load([load_address 'internal_polyhedron.mat'],'Rc_in')

    idx_eff=[];
    for jj=1:size(Rc,2)
        if inside_polyhedron_(Rc(:,jj),Rc_in)
            idx_eff=[idx_eff jj];
        end
    end
    idx_eff=idx_eff(~ismember(idx_eff,idx_wrong_disk));

    box_center=mean(rc_p,2);
    for jj=1:8
        [xi,yi,zi]=ind2sub([2,2,2],jj);
        idx_point=find((rc_p(1,:)-box_center(1))*(2*xi-3)>0&...
            (rc_p(2,:)-box_center(2))*(2*yi-3)>0&...
            (rc_p(3,:)-box_center(3))*(2*zi-3)>0);
        idx_particle=idx_contact(2,idx_point);
        idx_particle=intersect(idx_particle,idx_eff);

        Phi_tmp=length(idx_particle)*Vp/sum(Vcell(idx_particle));
        Phi_loop(jj,ii)=Phi_tmp;
        s2_tmp=mean(Ori(3,idx_particle).^2)*1.5-0.5;
        s2_loop(jj,ii)=s2_tmp;
        rsm_loop(jj,ii)=mean(rs(idx_point));

        idx_phi=find(Phi_tmp>=Phi_range,1,'last');
        rs_ALL{idx_phi}=[rs_ALL{idx_phi} rs(idx_point)];
        vcell_ALL{idx_phi}=[vcell_ALL{idx_phi} Vcell(idx_particle)];
    end
    disp(ii)
end



%% 
figure(1);clf;hold on
for ii=1:length(rs_ALL)
    [p,x]=hist(rs_ALL{ii},100);p=p/sum(p)/(x(2)-x(1));
    plot(x.^2,p)
end
legend


%%
Phi_loop=Phi_loop(:)';
rsm_loop=rsm_loop(:)';

Phi_range=linspace(min(Phi_loop),max(Phi_loop),10);
% Phi_range=linspace(min(Phi_loop),max(Phi_loop),7);
Phi_ave=zeros(1,length(Phi_range)-1);
Phi_std=zeros(1,length(Phi_range)-1);
rsm_ave=zeros(1,length(Phi_range)-1);
rsm_std=zeros(1,length(Phi_range)-1);
for ii=1:length(Phi_range)-1
    idx=find(Phi_loop>=Phi_range(ii)&Phi_loop<Phi_range(ii+1));
    Phi_ave(ii)=mean(Phi_loop(idx));
    Phi_std(ii)=std(Phi_loop(idx));
    rsm_ave(ii)=mean(rsm_loop(idx));
    rsm_std(ii)=std(rsm_loop(idx));
end


f=@(B,x)B*x.^(-1/3);
B=nlinfit(Phi_ave,rsm_ave,f,[1]);
figure(1);clf;hold on
plot(Phi_ave,rsm_ave,'ko')
xf=linspace(0.52,0.72,100);
yf=f(B,xf);
plot(xf,yf,'r-')


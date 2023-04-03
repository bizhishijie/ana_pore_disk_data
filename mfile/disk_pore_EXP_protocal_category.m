clear

exp_protocal={'a5','a5','a5','a5','a4','a4','a4','a4','a3','a3','a3','a3','a2','a2','a2',...
'a2','a1','a1','a1','a1','a5','a4','a3','a2','a1','a5','a4','a3','a2','a1',...
'c1','c2','c3','c4','c5','d1','d2','d3','d4','d5','c1','c2','c3','c4','c5',...
'b1','b2','b3','b4','b5','b1','b2','b1','b3','b4','b5','b3','b2','b4','b5','b2','b3',...
'b5','b2','b5','b3','b2','b3','b2','b3','b2','b2','b3','b5','b3'
};

D=30.75*2/0.8;
H=4.81*2/0.8;
Vp=D^2/4*H*pi;

figure(1);clf
hold on

Phi_all=[];
S2_all=[];

load_address=['D:\xiachj\research\src disk pore\data\'];
for ii=1:75
    load([load_address num2str(ii) '_disk\all_basic_data.mat'])
    Phi=length(idx_eff)*Vp/sum(Vcell(idx_eff));
    Phi_all=[Phi_all Phi];
    if Phi<0.5
%         disp(ii)
    end
    S2=mean(Ori(3,idx_eff).^2)*3/2-1/2;
    S2_all=[S2_all S2];
    exp_tmp=exp_protocal{ii};
    switch exp_tmp(1)
        case 'a'
            h=plot(Phi,S2,'ro');
            set(h,'MarkerSize',str2double(exp_tmp(2))*3)
        case 'b'
            h=plot(Phi,S2,'ks');
            set(h,'MarkerSize',str2double(exp_tmp(2))*3)
        case 'c'
            h=plot(Phi,S2,'b^');
        case 'd'
            h=plot(Phi,S2,'gv');
    end
end
axis([0.5 0.75 -0.2 1])

%%
idx_fixphi=find(Phi_all>0.65&Phi_all<0.68);
s2_range=linspace(min(S2_all(idx_fixphi)),max(S2_all(idx_fixphi))+eps,4);
for ii=1:length(s2_range)-1
    idx_tmp=S2_all(idx_fixphi)>=s2_range(ii)&S2_all(idx_fixphi)<s2_range(ii+1);
    idx_fixphi(idx_tmp)
end








% 读取圆盘位置
clear
% d=30.75*2/0.8;
% h=4.81*2/0.8;
% r=d/2;
% cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
load('../basic/pack_num_category.mat');
% load('p.mat');p=p';p_length=length(p);
hist_num=200;
for ii=1:length(pack_num_category)
    pack_num_category_tmp=pack_num_category{ii};
    x=linspace(0,8,hist_num);
    h=zeros(1,hist_num);
    for jj=1:length(pack_num_category_tmp)
        load(['..\basic\pack\pack' num2str(pack_num_category_tmp(jj)) '\rc_pore.mat'],'rc_pore');
        load(['..\basic\pack\pack' num2str(pack_num_category_tmp(jj)) '\pore_1.mat'],'is_del')
        h=h+hist(cellfun(@(c)mean(c),rc_pore(~is_del)),hist_num);
    end
    h=h/sum(h);
    plot(x,h,'-')
    hold on
end
legend('1','2','3','4');
set(gca,'yScale','log')
saveas(gcf,'..\fig\pore_rc.jpg')
saveas(gcf,'..\fig\pore_rc.fig')
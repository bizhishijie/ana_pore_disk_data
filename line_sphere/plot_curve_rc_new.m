% 读取圆盘位置
% 最终画的图已归一化
clear
fileList=dir('..\basic\pack\pack*');
load('../basic/pack_num_category.mat');
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%

clf
hist_num=100;
hist_length=80;
hist_centre=linspace(0,hist_length,hist_num)+hist_length/2/hist_num;
rc_num=zeros(4,length(pack_num_category));
for ii=1:length(pack_num_category)
    pack_num_category_tmp=pack_num_category{ii};
    %     clf
    rc_1_list=[];
    rc_2_list=[];
    rc_3_list=[];
    rc_a_list=[];
    for jj=1:length(pack_num_category_tmp)
        load(['..\basic\pack\pack' num2str(pack_num_category_tmp(jj)) '\rc_list.mat'],'rc_list')
        rc_1=rc_list(1,rc_list(2,:)==1);
        rc_2=rc_list(1,rc_list(2,:)==2);
        rc_3=rc_list(1,rc_list(2,:)==3);

        rc_num(1,ii)=rc_num(1,ii)+length(rc_1);
        rc_num(2,ii)=rc_num(2,ii)+length(rc_2);
        rc_num(3,ii)=rc_num(3,ii)+length(rc_3);
        rc_num(4,ii)=rc_num(4,ii)+size(rc_list,2);
        % 最小值取的三种情况

        rc_1_list=[rc_1_list rc_1];
        rc_2_list=[rc_2_list rc_2];
        rc_3_list=[rc_3_list rc_3];
        rc_a_list=[rc_a_list rc_list(1,:)];
    end

    [counts_1,~] = hist(rc_1_list,hist_centre);
    [counts_2,~] = hist(rc_2_list,hist_centre);
    [counts_3,~] = hist(rc_3_list,hist_centre);
    [counts_a,~] = hist(rc_a_list,hist_centre);

    figure(1)
    hold on
    plot(hist_centre,counts_1/rc_num(1,ii));
    %     plot(hist_centre.^2,counts_1/rc_num(1,ii));
    figure(2)
    hold on
    plot(hist_centre,counts_2/rc_num(2,ii));
    %     plot(hist_centre.^2,counts_2/rc_num(2,ii));
    figure(3)
    hold on
    plot(hist_centre,counts_3/rc_num(3,ii));
    %     plot(hist_centre.^2,counts_3/rc_num(3,ii));
    figure(4)
    hold on
    plot(hist_centre,counts_a/rc_num(4,ii));
    %     plot(hist_centre.^2,counts_a/rc_num(4,ii));
end
figure(1)
set(gca,'yScale','log')
% legend('角','底','面','总')
legend('1','2','3','4')
title('角')
box on
saveas(gcf,'..\fig\rc_pack_角.jpg')
saveas(gcf,'..\fig\rc_pack_角.fig')

figure(2)
set(gca,'yScale','log')
% legend('角','底','面','总')
legend('1','2','3','4')
title('底')
box on
saveas(gcf,'..\fig\rc_pack_底.jpg')
saveas(gcf,'..\fig\rc_pack_底.fig')

figure(3)
set(gca,'yScale','log')
% legend('角','底','面','总')
legend('1','2','3','4')
title('面')
box on
saveas(gcf,'..\fig\rc_pack_面.jpg')
saveas(gcf,'..\fig\rc_pack_面.fig')

figure(4)
set(gca,'yScale','log')
% legend('角','底','面','总')
legend('1','2','3','4')
title('总')
box on
saveas(gcf,'..\fig\rc_pack_总.jpg')
saveas(gcf,'..\fig\rc_pack_总.fig')
% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
load('../basic/pack_num_category.mat');
% d=30.75*2;
% h=4.81*2;
% r=d/2;
% cylinder_size=[r,h];

clf
hist_num=30;
hist_centre=linspace(0,70,hist_num);
for ii=1:length(pack_num_category)
    pack_num_category_tmp=pack_num_category{ii};
    %     clf
    counts_1_list=zeros(hist_num,length(pack_num_category_tmp));
    counts_2_list=zeros(hist_num,length(pack_num_category_tmp));
    counts_3_list=zeros(hist_num,length(pack_num_category_tmp));
    counts_a_list=zeros(hist_num,length(pack_num_category_tmp));
    for jj=1:length(pack_num_category_tmp)
        load(['..\basic\pack\pack' num2str(pack_num_category_tmp(jj)) '\rc_list.mat'],'rc_list')
        rc_list_1=rc_list(:,rc_list(2,:)==1);
        rc_list_2=rc_list(:,rc_list(2,:)==2);
        rc_list_3=rc_list(:,rc_list(2,:)==3);

        [counts_1,~] = hist(rc_list_1(1,:),hist_centre);
        [counts_2,~] = hist(rc_list_2(1,:),hist_centre);
        [counts_3,~] = hist(rc_list_3(1,:),hist_centre);
        [counts_a,~] = hist(rc_list(1,:),hist_centre);

        %         plot(centers_1,counts_1,'r')
        %         plot(centers_2,counts_2,'g')
        %         plot(centers_3,counts_3,'b')
        %         plot(centers_a,counts_a,'k')

        counts_1_list(:,jj)=counts_1;
        counts_2_list(:,jj)=counts_2;
        counts_3_list(:,jj)=counts_3;
        counts_a_list(:,jj)=counts_a;
    end
    figure(1)
    hold on
    %     plot(hist_centre.^2,mean(counts_1_list,2));
    plot(hist_centre,mean(counts_1_list,2));
    figure(2)
    hold on
    %     plot(hist_centre.^2,mean(counts_2_list,2));
    plot(hist_centre,mean(counts_2_list,2));
    figure(3)
    hold on
    %     plot(hist_centre.^2,mean(counts_3_list,2));
    plot(hist_centre,mean(counts_3_list,2));
    figure(4)
    hold on
    %     plot(hist_centre.^2,mean(counts_a_list,2));
    plot(hist_centre,mean(counts_a_list,2));
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
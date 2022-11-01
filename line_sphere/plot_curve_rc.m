% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
d=30.75*2;
h=4.81*2;
r=d/2;
cylinder_size=[r,h];

for ii=1:length(fileList)
    load(['..\basic\pack\' fileList(ii).name '\rc_list.mat'],'rc_list')
    clf
    hold on
    rc_list_1=rc_list(:,rc_list(2,:)==1);
    rc_list_2=rc_list(:,rc_list(2,:)==2);
    rc_list_3=rc_list(:,rc_list(2,:)==3);
    [counts_1,centers_1] = hist(rc_list_1(1,:),50);
    [counts_2,centers_2] = hist(rc_list_2(1,:),50);
    [counts_3,centers_3] = hist(rc_list_3(1,:),50);
    [counts_a,centers_a] = hist(rc_list(1,:),50);
    plot(centers_1,counts_1)
    plot(centers_2,counts_2)
    plot(centers_3,counts_3)
    plot(centers_a,counts_a)
    set(gca,'yScale','log')
    legend('角','底','面','总')
    box on
    title(fileList(ii).name)
    saveas(gcf,['..\fig\rc_pack_' fileList(ii).name '.jpg'])
    %     save(['..\basic\pack\' fileList(ii).name '\rc_list_hist.mat'],'centers_1','counts_2','centers_3','centers_a' ...
    %         ,'counts_1','counts_2','counts_3','counts_a');
end
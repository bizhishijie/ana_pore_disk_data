% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
for ii=1:length(fileList)
    clf
    load(['..\basic\pack\' fileList(ii).name '\length.mat'],'length_list','line_forward_list')
    [counts,centers]=hist(length_list,100);
    plot(centers,counts)
    set(gca,'yScale','log')
    box on
    title(fileList(ii).name)
    saveas(gcf,['..\fig\line_pack_' fileList(ii).name '.jpg'])
end
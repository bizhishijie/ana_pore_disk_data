% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
for ii=1:length(fileList)
    clf
    load(['..\basic\pack\' fileList(ii).name '\length.mat'],'length_list','line_forward_list')

    angle_num=100;
    angle_list=linspace(-1,1-2/angle_num,angle_num);
    hist_angle=zeros(1,angle_num);
    for jj =1:length(angle_list)
        tmp=length_list.*(line_forward_list(3,:)>=angle_list(jj)& ...
            line_forward_list(3,:)<(angle_list(jj)+2/angle_num));
        hist_angle(jj)=mean(tmp(tmp~=0));
    end

    polarplot(acos(angle_list),hist_angle);
    title(fileList(ii).name)
    saveas(gcf,['..\fig\line_angle_pack_' fileList(ii).name '.jpg'])
end
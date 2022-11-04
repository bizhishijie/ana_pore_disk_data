
clear
fileList=dir('..\basic\pack\pack*');
load('pack_num_category.mat');
angle_num=150;
angle_list=linspace(0,1,angle_num);
for ii=1:length(pack_num_category)
    pack_num_category_tmp=pack_num_category{ii};
    %     clf
    %     figure(ii)
    hist_angle=zeros(length(pack_num_category_tmp),angle_num);
    for jj=1:length(pack_num_category_tmp)
        load(['..\basic\pack\pack' num2str(pack_num_category_tmp(jj)) '\length.mat'],'length_list','line_forward_list')
        for kk =1:length(angle_list)
            tmp=length_list.*(abs(line_forward_list(3,:))>=(angle_list(kk)-1/2/angle_num)& ...
                abs(line_forward_list(3,:))<(angle_list(kk)+1/2/angle_num));
            hist_angle(jj,kk)=mean(tmp(tmp~=0));
        end
    end
    polarplot(acos(angle_list),mean(hist_angle,1));
    hold on
    %     title(num2str(ii))
end
legend('1','2','3','4')
axis([0 90 0 60])
box on
saveas(gcf,'..\fig\line_pack_angle.jpg')
saveas(gcf,'..\fig\line_pack_angle.fig')
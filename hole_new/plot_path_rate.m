clear
fileList=dir('..\basic\pack\pack*');
load('../basic/pack_num_category.mat');
angle_num=90;
angle_list=linspace(0,1,angle_num);
for ii=1:length(pack_num_category)
    pack_num_category_tmp=pack_num_category{ii};
    %     clf
    %     figure(ii)
    hist_angle=zeros(length(pack_num_category_tmp),angle_num);
    for jj=1:length(pack_num_category_tmp)
        rate_list=load(['..\basic\pack\pack' num2str(pack_num_category_tmp(jj)) '\path_rate.mat'],'rate_list').rate_list;
        line_forward_list=load(['..\basic\pack\pack' num2str(pack_num_category_tmp(jj)) '\path_rate.mat'],'p_list').p_list;
        line_forward_list=line_forward_list(:,rate_list>1);
        rate_list=rate_list(rate_list>1);
        line_forward_list=(line_forward_list(1:3,:)-line_forward_list(4:6,:));
        line_forward_list=line_forward_list./vecnorm(line_forward_list);
        for kk =1:length(angle_list)
            tmp=rate_list.*(abs(line_forward_list(3,:))>=(angle_list(kk)-1/2/angle_num)& ...
                abs(line_forward_list(3,:))<(angle_list(kk)+1/2/angle_num));
            hist_angle(jj,kk)=mean(tmp(tmp~=0));
        end
    end
    polarplot(acos(angle_list),mean(hist_angle,1));
    hold on
    %     title(num2str(ii))
end
legend('1','2','3','4')
axis([0 90 0 3])
box on
saveas(gcf,'..\fig\path_rate.jpg')
saveas(gcf,'..\fig\path_rate.fig')
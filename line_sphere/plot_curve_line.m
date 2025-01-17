% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
load('../basic/pack_num_category.mat');
clf
hold on

hist_num=500;max_lengh=80;
hist_centre=linspace(0+max_lengh/hist_num/2,max_lengh-max_lengh/hist_num/2,hist_num);
for ii=1:length(pack_num_category)
    pack_num_category_tmp=pack_num_category{ii};
    length_tmp_list=zeros(hist_num,length(pack_num_category_tmp));
    for jj=1:length(pack_num_category_tmp)
        load(['..\basic\pack\pack' num2str(pack_num_category_tmp(jj)) '\length.mat'],'length_list','line_forward_list')
        [counts,~]=hist(length_list,hist_centre);
        %         plot(hist_centre,counts)
        %         set(gca,'yScale','log')
        %         box on
        %         title(fileList(ii).name)
        length_tmp_list(:,jj)=counts;
    end
    tmp=mean(length_tmp_list,2);
    tmp(end)=0;
    plot(hist_centre,tmp/sum(tmp));
end

set(gca,'yScale','log')
legend('1','2','3','4')
box on
saveas(gcf,'..\fig\line_pack.jpg')
saveas(gcf,'..\fig\line_pack.fig')
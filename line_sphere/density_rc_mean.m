clear;clc
load('pack_density.mat')
load('pack_num_category.mat')
fileList=dir('..\basic\pack\pack*');
file_num=zeros(1,length(fileList));
rc_mean=cell(1,length(pack_num_category));

for ii=1:length(fileList)
    file_num(ii)=str2num(fileList(ii).name(5:end));
end
file_num=sort(file_num)';
for ii=1:length(pack_num_category)
    pack_list_tmp=pack_num_category{ii};
    rc_mean{ii}=zeros(1,length(pack_list_tmp));
    for jj=1:length(pack_list_tmp)
        load(['..\basic\pack\pack' num2str(pack_list_tmp(jj)) '\rc_list.mat'],'rc_list')
        rc_mean{ii}(jj)=mean(rc_list(1,:));
    end
end
hold on
for ii=1:length(pack_num_category)
    plot(pack_density( ismember(file_num,pack_num_category{ii}) ).^(1/3),1./rc_mean{ii},'o')
end
legend('1','2','3','4')

% plot(pack_density.^(1/3),1./rc_mean,'o')
plot([0 ,mean(pack_density.^(1/3))],[0 ,mean(1./cell2mat(rc_mean))]);

load('pack_density.mat')
load("file_name.mat")
load('..\basic\pack_num_category.mat')
[~,min_id]=min(pack_density);
[~,max_id]=max(pack_density);
hist_num=100;
hist_length=80;
hist_centre=linspace(0,hist_length,hist_num)+hist_length/2/hist_num;
rc_num=zeros(4,length(pack_num_category));

fileList=dir('..\basic\pack\pack*');

load(['..\basic\pack\pack' num2str(file_name(min_id)) '\rc_list.mat'],'rc_list')
rc_1=rc_list(1,rc_list(2,:)==1);
rc_2=rc_list(1,rc_list(2,:)==2);
rc_3=rc_list(1,rc_list(2,:)==3);
rc_a=rc_list(1,:);
% 最小值取的三种情况

[counts_1,~] = hist(rc_1,hist_centre);
[counts_2,~] = hist(rc_2,hist_centre);
[counts_3,~] = hist(rc_3,hist_centre);
[counts_a,~] = hist(rc_a,hist_centre);
figure(1)
clf
hold on
% plot(hist_centre,counts_1,'ro-');
% plot(hist_centre,counts_2,'go-');
% plot(hist_centre,counts_3,'bo-');
plot(hist_centre.^2,counts_a,'ko-');

set(gca,'yScale','log')
%%
load(['..\basic\pack\pack' num2str(file_name(max_id)) '\rc_list.mat'],'rc_list')
rc_1=rc_list(1,rc_list(2,:)==1);
rc_2=rc_list(1,rc_list(2,:)==2);
rc_3=rc_list(1,rc_list(2,:)==3);
rc_a=rc_list(1,:);
% 最小值取的三种情况

[counts_1,~] = hist(rc_1,hist_centre);
[counts_2,~] = hist(rc_2,hist_centre);
[counts_3,~] = hist(rc_3,hist_centre);
[counts_a,~] = hist(rc_a,hist_centre);
% figure(2)
% clf
% hold on
% plot(hist_centre,counts_1,'r*-');
% plot(hist_centre,counts_2,'g*-');
% plot(hist_centre,counts_3,'b*-');
plot(hist_centre.^2,counts_a,'k*-');

set(gca,'yScale','log')
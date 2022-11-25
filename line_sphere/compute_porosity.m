% 读取圆盘位置
clear
fileList=dir('..\basic\pack\pack*');
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
% cylinder_size=[r,h];%
dz=10;% 统计圆盘质心位于高度在z~z+dz的孔隙度的占比
epsilon_list=zeros(1,length(fileList));
S_list=zeros(1,length(fileList));
for ii=1:length(fileList)
    Rc=load(['..\basic\pack\' fileList(ii).name '\basic.mat']).Rc;
    Ori=load(['..\basic\pack\' fileList(ii).name '\basic.mat']).Ori;
    edge=load(['..\basic\pack\' fileList(ii).name '\edge.mat']).edge;
    is_in=load(['..\basic\pack\' fileList(ii).name '\rc_is_in.mat']).is_in;

    disp(fileList(ii).name)

    Rc_inner=Rc(:,is_in);
    Ori_inner=Ori(:,is_in);
    %     show_cylinder(Rc_inner,Ori_inner);

    Rc_z_max=max(Rc_inner(3,:));
    Rc_z_min=min(Rc_inner(3,:));

    z_array=Rc_z_min:dz:Rc_z_max;
    epsilon=zeros(size(z_array));
    n=zeros(size(z_array));
    for jj=1:length(z_array)
        z=z_array(jj);
        Rc_tmp=Rc(:,(Rc(3,:)<z+dz)&(Rc(3,:)>=z));
        n(jj)=length(Rc_tmp);

        Rc_max=max(Rc_tmp,[],2);
        Rc_min=min(Rc_tmp,[],2);
        volume_all=(Rc_max(2)-Rc_min(2))*(Rc_max(1)-Rc_min(1))*dz;
        cyliner_num=sum(any(Rc_tmp));
        epsilon(jj)=cyliner_num*pi*r^2*h/volume_all;
    end

    %     plot(epsilon,z_array)
    %     axis([0 1 min(z_array) max(z_array)])

    epsilon_list(ii)=mean(epsilon(floor(end/4):ceil(end*3/4)));
    S=mean(3*sum( Ori_inner .* repmat([0;0;1],1,size(Ori_inner,2) ) ).^2 - 1 )/2;
    S_list(ii)=S;
end

load('pack_num_category.mat')
file_num=zeros(1,length(fileList));
for ii=1:length(fileList)
    file_num(ii)=str2num(fileList(ii).name(5:end));
end
file_num=sort(file_num)';
hold on
for ii =1:length(pack_num_category)
    file_id=ismember(file_num,pack_num_category{ii});
    plot(S_list(file_id),epsilon_list(file_id),'o')
end
legend('1','2','3','4')
xlabel('S')
ylabel('\epsilon')
axis([0 1 0 1])
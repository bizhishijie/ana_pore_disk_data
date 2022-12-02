clear
d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
load('p.mat')
p=p';
for nn=1:length(fileList)
    load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
    p_all=cell(length(Rc),1);
    normal_vec=[0;0;1];
    length_p=size(p,2);
    length_Rc=size(Rc,2);
    for ii=1:size(Rc,2)
        % 已知变换前后的法向量，求得变换矩阵
        v=cross(normal_vec,Ori(:,ii));
        %     s=norm(v);
        c=dot(normal_vec,Ori(:,ii));
        v=[0        -v(3)   v(2);
            v(3)     0       -v(1);
            -v(2)    v(1)    0];
        R=diag(ones(1,3))+v+v*v/(1+c);
        p_new=R*p+Rc(:,ii);
        p_all{ii}=p_new;

        %         d=30.75*2/0.8;
        %         h=4.81*2/0.8;
        %         r=d/2;
        %         cylinder_size=[r,h];

        %         draw_cylinder(Rc(:,ii),Ori(:,ii),'r',r,h,0.5)
        %         plot3(p_new(1,:), ...
        %             p_new(2,:), ...
        %             p_new(3,:),'.')
        % %         验证旋转的正确
    end
    save(['..\basic\pack\' fileList(nn).name '\p_all.mat'],'p_all')
end
function [V,R]=cylinder_to_voronoi(Rc,Ori,p)
% Rc是质心坐标
% Ori是朝向
% 返回的是晶胞的顶点，区域，和有效点的id
% tri 三角划分是 distmesh 文件夹下的 tri_cylinder 文件划分的
% 也可以根据 cylinder_sample_generate 划分圆柱的采样点
% 旋转算法见 https://zhuanlan.zhihu.com/p/134299994

normal_vec=[0;0;1];
length_p=size(p,2);
length_Rc=size(Rc,2);
p_all=zeros(3,length_p*length_Rc);

for ii=1:size(Rc,2)
    % 已知变换前后的法向量，求得变换矩阵
    v=cross(normal_vec,Ori(:,ii));
    %     s=norm(v);
    c=dot(normal_vec,Ori(:,ii));
    v=[0        -v(3)   v(2);
        v(3)     0       -v(1);
        -v(2)    v(1)    0];
    R=diag(ones(1,3))+v+v*v/(1+c);
    p_all(:,(ii-1)*length_p+1:ii*length_p)=R*p+Rc(:,ii);

%     d=30.75*2/0.8;
%     h=4.81*2/0.8;
%     r=d/2;
%     cylinder_size=[r,h];
% 
%     draw_cylinder(Rc(:,ii),Ori(:,ii),'r',r,h,0.5)
%     plot3(p_all(1,(ii-1)*length_p+1:ii*length_p), ...
%         p_all(2,(ii-1)*length_p+1:ii*length_p), ...
%         p_all(3,(ii-1)*length_p+1:ii*length_p),'.')


    %         验证旋转的正确

end

dt = delaunayTriangulation(p_all');
[V,R] = voronoiDiagram(dt);
end
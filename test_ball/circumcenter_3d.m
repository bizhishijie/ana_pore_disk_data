function [circumcenter] = circumcenter_3d(p1, p2, p3)
% 计算三维空间中三角形的外心
% 输入：三个点的坐标 p1, p2, p3 (3x1向量)
% 输出：外心坐标 circumcenter (3x1向量)

% 计算三条边的中点
p12 = (p1 + p2) / 2;
p13 = (p2 + p3) / 2;

% 计算三条边的向量
v12 = p2 - p1; v12=v12/norm(v12);
v13 = p3 - p2; v13=v13/norm(v13);

n = cross(v12,v13);n=n/norm(n);
n12 = cross(v12,n);
n13 = cross(v13,n);

k=[1,-dot(n12,n13);dot(n12,n13),-1]\[dot(p13-p12,n12);dot(p13-p12,n13)];
circumcenter=k(1)*n12+p12;
end
function cross_p=intersects_line_flat(line_p,line_forward,flat_p,flat_normal)
% 返回直线和平面的交点
% 参数分别为直线上的点，直线的方向向量，平面上的点，平面的法向量
% https://wuli.wiki/online/LPint.html
% intersects_line_flat([1;1;1],[1;1;1],[1;1;1],[1;1;1])

cross_p=dot((flat_p-line_p),flat_normal)/dot(line_forward,flat_normal)*line_forward+line_p;
end
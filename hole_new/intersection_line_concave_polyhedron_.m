function p_inter=intersection_line_concave_polyhedron_(p,f,v,k)
% 射线上点p,射线的方向f,多面体的点的集合v,点的连接方式k
% 计算射线和凹多面体的交点
% intersection_line_concave_polyhedron_([100,100,100],[0,0,1],v,ks{1})
%%
p_inter=[];
p_max=max(v,2);
for ii=1:length(k)
    p_inter=[p_inter intersection_line_face_([p;p+p_max.*f]',v(k(ii,:),:))];
end
end
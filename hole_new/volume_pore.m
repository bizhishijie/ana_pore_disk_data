function v_pore=volume_pore(ks_p,pore_c_id_t,v,Rc,Ori,p)
Rc_t=Rc(:,pore_c_id_t);
Ori_t=Ori(:,pore_c_id_t);
%%
% tri=triangulation(ks_p-1,v(2:end,:));% 将其与co求交集，计算体积后相减
% 输入的是面的tri基本没啥用
% 避开inf的问题
% tetra=cell(1,length(ks_p));
% parfor ii=1:length(ks_p)
%     new_edge=setdiff(ks_p(sum(ismember(ks_p,ks_p(ii,:)),2)==2,:),ks_p(ii,:));
%     tetra{ii}=[repmat(ks_p(ii,:),size(new_edge,1),1) new_edge];
% end
% tetra=unique(sort(cell2mat(tetra'),2),'rows');
% v_all=zeros(1,length(tetra));
% v_inter=zeros(1,length(tetra));
% parfor ii=1:length(tetra)
%     warning('off');
%     try
%         [k,v_all(ii)]=convhulln(v(tetra(ii,:),:));
%     catch
%         [k,v_all(ii)]=convhulln(v(tetra(ii,:),:),{'QJ'});
%     end
%     for jj=1:size(Rc_t,2)
%         v_inter(ii)=v_inter(ii)+intersection_polyhedron2_(k,pco{jj},v(tetra(ii,:),:)',pc{jj});
%     end
% end
% 太慢了，时间不可接受
%%
P=v(unique(ks_p),:);
k = boundary(P);
P=P(unique(k(:)),:);
try
    [shp,v_all] = convhulln(P);% 凹的部分和颗粒接触，忽略吧
catch
    v_pore=0;
    return
end
v_inter=0;
for ii=1:size(Rc_t,2)
    Rc_p=cylinder_to_point(Rc_t(:,ii),Ori_t(:,ii),p);
    v_inter=v_inter+intersection_polyhedron2_(shp,convhulln(Rc_p'),P',Rc_p);
end
v_pore=v_all-v_inter;
if v_pore<0
    v_pore=0;
end
end
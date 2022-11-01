for jj=1:size(Rc,2)
    % 对每一个圆盘遍历，对应regions中的前若干个点
    id_tmp=regions((jj-1)*p_length+1:jj*p_length);
    id_tmp=cell2mat(cellfun(@transpose,id_tmp,'UniformOutput',false));
    id_tmp=sort(id_tmp);
    id_uni=unique(id_tmp);
    id_statistics=zeros(length(id_uni),2);
    id_statistics(:,1)=id_uni;
    for kk=1:size(id_uni)
        id_statistics(kk,2)=sum(id_tmp==id_uni(kk));
    end
    id_uni(id_statistics(:,2)>=4)=[];
    % 在 id  _tmp 里面，出现4次及以上的是内部的点，应该删除
    cylinder_point=vertice(id_uni,:);
    % 找到对应的点的编号
    try
        K = convhull(cylinder_point);
    catch
        continue
    end
    defaultFaceColor  = [0.6875 0.8750 0.8984];
    trisurf(K, cylinder_point(:,1),cylinder_point(:,2),cylinder_point(:,3) , ...
        'FaceColor', defaultFaceColor, 'FaceAlpha',0.8)
    hold on
    disp(jj)
end
show_cylinder(Rc,Ori)
axis([0 500 0 500 0 500])
saveas(gcf,'test.jpg')
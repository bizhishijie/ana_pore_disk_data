%% 计算每一个孔的path到周围三个圆柱的最近点，放到列表中
% for ii=1:length(pore)
%     pore_edge=[];
%     if pore(ii).is_del
%         continue
%     end
%     path_cell=pore(ii).path;
%     for jj=1:length(path_cell)% 对path遍历
%         % 对应的就是第jj个邻居
%         common_cylinder=intersect(pore(pore(ii).neighbor(jj)).near_cylinder,pore(ii).near_cylinder);
%         % 当邻居数目较多的时候，共同的圆柱不一定只有三个
%         path_tmp=path_cell{jj};
%         path_tmp(4,:)=sqrt(sum((path_tmp(1:3,:)-pore(ii).rc).^2));% path的第四行是该点到rc的距离，根据它排序
%         path_tmp=sortrows(path_tmp',4)';
%         path_tmp=path_tmp(1:3,:);% 删除该点到rc的距离这一行
%
%         dis=inf(1,length(path_tmp));
%         for kk=1:size(path_tmp,2)
%             if ~is_point_in_cylinder(path_tmp(:,kk),Rc,Ori,cylinder_size)
%                 dis_tmp=distance_point_cylinder(path_tmp(:,kk),Rc,Ori,cylinder_size);
%                 dis(kk)=min(dis_tmp(1,:));
%             end
%         end
%         [~,idx]=min(dis);
%         pore_edge=[];
%         for kk=1:idx% 对path上的点遍历
%             for cc=common_cylinder'
%                 Ac=[Rc(:,cc)-Ori(:,cc)*h/2 Rc(:,cc)+Ori(:,cc)*h/2];
%                 [Distance,P0,P1] = distance_cylinder_sphere(Ac, r, path_tmp(:,kk), 0);
%                 pore_edge=[pore_edge P0];
%             end
%             % 根据多个点拟合平面
%         end
%     end
%     pore(ii).edge=pore_edge;
% end
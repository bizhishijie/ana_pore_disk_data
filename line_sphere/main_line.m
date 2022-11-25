% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
line_num=1e5;
for ii=1:length(fileList)
    length_list=[];
    line_forward_list=[];
    load(['..\basic\pack\' fileList(ii).name '\basic.mat'])
    load(['..\basic\pack\' fileList(ii).name '\rc_is_in.mat'])
    disp(fileList(ii).name)

    Rc_max=max(Rc,[],2)-r;
    Rc_min=min(Rc,[],2)+r;
    Rc=Rc(:,is_in);
    %     Rc_outside=find(any((Rc>Rc_max)|(Rc<Rc_min)));
    %     Rc(:,Rc_outside)=[];
    %     Ori(:,Rc_outside)=[];% 去除比较靠外的圆盘

    parfor line_id=1:line_num% 改变循环次数
        %         Rc_cross_id_list=[];
        cross_p_list=[];
        cross_p_close_id=[];
        while size(cross_p_list,2)<4   % 交点太少了就重新随机，只和一个圆盘相交算不了
            line_p=rand(3,1).*(Rc_max-Rc_min)+Rc_min;
            line_forward=random_unit_vector;

            for jj=1:length(Rc)
                Rc_tmp=Rc(:,jj);
                Ori_tmp=Ori(:,jj);
                cross_p=intersects_line_cylinder(line_p,line_forward,Rc_tmp,Ori_tmp,cylinder_size);% 获得直线与圆柱体的交点
                if ~isempty(cross_p)
                    cross_p_list=[cross_p_list [cross_p;Ori_tmp,Ori_tmp]];% 交点的列表,以及对应的圆盘的法向量
                    %                 Rc_cross_id_list=[Rc_cross_id_list jj];
                end
            end
        end
        cross_p_list=sortrows(cross_p_list')';% 必然为偶数个交点，将其按照空间位置排序

        % 处理穿模
        cross_p_through_id=find(sum(cross_p_list(4:6,3:2:end)-cross_p_list(4:6,2:2:end-2),1)==0);% 找到穿模的点，根据法向量判断是否穿模
        cross_p_list(:,[2*cross_p_through_id 2*cross_p_through_id+1])=[];

        % 处理贴的过近的
        cross_p_close_id=find(abs(sum(cross_p_list(4:6,3:2:end).*cross_p_list(4:6,2:2:end-2)))>0.995&...
            dot(cross_p_list(1:3,3:2:end)-cross_p_list(1:3,2:2:end-2),cross_p_list(4:6,3:2:end))<h*0.05);

        %         if ~isempty(cross_p_close_id)
        %             disp(length(cross_p_close_id))
        %             show_cylinder(Rc(:,Rc_cross_id_list),Ori(:,Rc_cross_id_list));
        %             hold on
        %             p_1=line_p+1000*line_forward;p_2=line_p-1000*line_forward;
        %             line([p_1(1) p_2(1)],[p_1(2) p_2(2)],[p_1(3) p_2(3)])
        %             plot3(cross_p_list(1,2*cross_p_close_id),cross_p_list(2,2*cross_p_close_id),cross_p_list(3,2*cross_p_close_id),'b*');
        %             plot3(cross_p_list(1,2*cross_p_close_id+1),cross_p_list(2,2*cross_p_close_id+1),cross_p_list(3,2*cross_p_close_id+1),'b*');
        %         end
        % %         运行这一段需要Rc_cross_id_list

        cross_p_list(:,[2*cross_p_close_id 2*cross_p_close_id+1])=[];% 如果连续的

        cross_dis_tmp=cross_p_list(1:3,3:2:end)-cross_p_list(1:3,2:2:end-2);% 其中第偶数个为空隙的向量
        cross_dis_tmp=sqrt(sum(cross_dis_tmp.^2));
        for jj=1:(size(cross_p_list,2)/2)
            length_list=[length_list cross_dis_tmp];
            line_forward_list=[line_forward_list repmat(line_forward,1,size(cross_dis_tmp,2))];
        end

        %         n_length=200;% 画的直线的长度
        %         draw_cylinder(Rc_tmp,Ori_tmp,'r',r,h,1)
        %         p_1=line_p-n_length*line_forward;
        %         p_2=line_p+n_length*line_forward;
        %         line([p_1(1) p_2(1)],[p_1(2) p_2(2)],[p_1(3) p_2(3)])
        %         axis equal
        %         for ii=1:size(cross_p,2)
        %             plot3(cross_p(1,ii),cross_p(2,ii),cross_p(3,ii),'b*');
        %         end
    end
    save(['..\basic\pack\' fileList(ii).name '\length.mat'],'length_list','line_forward_list')
end

% histogram(length_list)
% set(gca,'yScale','log')

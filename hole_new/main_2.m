% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
load('p.mat');p=p';p_length=length(p);
for nn=1:15
    %% 加载数据
    disp(nn)
    load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
    load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
    load(['..\basic\pack\' fileList(nn).name '\edge.mat'])
    load(['..\basic\pack\' fileList(nn).name '\voro.mat'])
    %% 计算孔的属性
    pore_c_id=zeros(size(pore_rc,1),4);% 每个孔相邻的圆柱
    tbl=tabulate(cell2mat(v_cell'));
    parfor pp=1:size(pore_rc,1)
        pore_c_id(pp,:)=find(cellfun(@(v)ismember(map(pp),v),v_cell));
    end

    neighbor_cell=cell(1,size(pore_rc,1));% 孔的邻居
    parfor pp=1:size(pore_rc,1)
        neighbor_cell{pp}=find(sum(ismember(pore_c_id,pore_c_id(pp,:)),2)==3);
    end
    %%
    path_cell=cell(1,size(pore_rc,1));% 孔的路径
    throat_cell=cell(1,size(pore_rc,1));% 嗓子的位置
    parfor pp=1:size(pore_rc,1)
        path_cell{pp}=cell(1,length(neighbor_cell{pp}));
        throat_cell{pp}=cell(1,length(neighbor_cell{pp}));
        for jj=1:length(neighbor_cell{pp})% 对每个孔和它的邻居遍历
            qq=neighbor_cell{pp}(jj);
            cc=intersect(pore_c_id(pp,:),pore_c_id(qq,:));% common cylinder
            tbl=tabulate(cell2mat(v_cell(cc)'));
            path_tmp=tbl(tbl(:,2)==3,1);% 需要根据远近将path排序，便于合并时连接

            dis=sqrt(sum((v(path_tmp,:)-pore_rc(pp,:)).^2,2));
            path_tmp=sortrows([dis,path_tmp]);
            path_tmp=path_tmp(:,2);
            path_cell{pp}{jj}=path_tmp;% 排序后放入cell

            Rc_t=Rc(:,cc);Ori_t=Ori(:,cc);
            dl=zeros(1,length(path_tmp));% 根据距离最短找到嗓子的位置
            for ii=1:length(path_tmp)
                pa=path_tmp(ii);
                dis=distance_point_cylinder(v(pa,:),Rc_t',Ori_t',cylinder_size);
                dis=min(dis(1,:));
                dl(ii)=dis;
            end
            [~,idx]=min(dl);
            dp=diff(v(path_tmp,:));
            if idx>size(dp,1)% 如果第一个点就是throat，实际上可能是由于采样点数目不够多造成的
                dp=dp(idx-1,:);
            else
                dp=dp(idx,:);
            end
            throat_cell{pp}{jj}=[path_tmp(idx) dp/norm(dp)];
            % 第一个数是嗓子所在点的编号，后三个数是path在嗓子处的朝向，朝向总是向外的
        end
    end
    %% 计算孔拥有的采样点
    ks=cell(1,size(pore_rc,1));% 孔有的采样点和连接
    parfor pp=1:size(pore_rc,1)
        ks{pp}=[];
        for jj=1:4
            ctmp=pore_c_id(pp,jj);% cylinder
            stmp=find(cell2mat(cellfun(@(c)any(ismember(c,pore_c_id(pp,:))),co{ctmp},'UniformOutput',false)));% 孔包含的采样点
            sptmp=k_cell{ctmp}(stmp);
            for ii=1:length(stmp)
                p_id=unique(sptmp{ii});
                sutmp=v(p_id,:);% 获得了边界点，接下来需要根据throat的位置和方向删除
                flag=true(size(sutmp,1),1);
                for tt=1:length(throat_cell{pp})
                    s_2=side_of_flat(sutmp,v(throat_cell{pp}{tt}(1),:),throat_cell{pp}{tt}(2:4));
                    s_1=side_of_flat(pore_rc(pp,:),v(throat_cell{pp}{tt}(1),:),throat_cell{pp}{tt}(2:4));
                    if s_1==0
                        continue
                    end
                    flag=flag & (s_1==s_2);
                end
                p_id=p_id(flag,:);% 找到了在内侧的点的id
                ks{pp}=[ks{pp};sptmp{ii}(all(ismember(k_cell{ctmp}{stmp(ii)},p_id),2),:)];
            end
        end
    end
    %% 展示pore
    %     clf
    %     pp=120;
    %     show_pore(pore_rc(pp,:),ks{pp},v,Rc(:,pore_c_id(pp,:)),Ori(:,pore_c_id(pp,:)),pore_rc(neighbor_cell{pp},:),path_cell(pp),throat_cell{pp});
    %% 根据文章定义rk并合并
    % pore的属性有rk path_cell ks neighbor_cell
    %     save
    flag=true;
    is_del=false(1,size(pore_rc,1));
    cnt=1;
    pore_c_id=num2cell(pore_c_id,2);
    while flag
        flag=false;
        disp(cnt)
        cnt=cnt+1;
        for pp=1:size(pore_rc,1)
            if ~is_del(pp)% 如果没删去
                for jj=1:length(neighbor_cell{pp})
                    qq=neighbor_cell{pp}(jj);
                    if ~is_del(qq)% 如果没删去
                        if norm(pore_rc(pp,:)-pore_rc(qq,:))<rk(pp) && rk(qq)<rk(pp)
                            % 满足条件合并
                            is_del(qq)=true;
                            flag=true;

                            nei_t1=neighbor_cell{pp};
                            idx_p=find(nei_t1==qq);
                            nei_t1(nei_t1==qq)=[];
                            nei_t2=neighbor_cell{qq};
                            idx_q=find(nei_t2==pp);
                            nei_t2(nei_t2==pp)=[];
                            an=setdiff(nei_t2,nei_t1);% additional neighbor
                            an(is_del(an))=[];% 已经被删除的邻居就不算了

                            if ~isempty(an)% 如果有新的邻居
                                path_t=path_cell{pp}{idx_p};% 将要被合并的path
                                path_cell{pp}(idx_p)=[];
                                neighbor_cell{pp}(idx_p)=[];
                                for kk=1:length(neighbor_cell{qq})
                                    nk=neighbor_cell{qq}(kk);
                                    if ismember(nk,an)
                                        % 添加路径
                                        path_cell{pp}{end+1}=[path_t(1:end-1);path_cell{qq}{kk}];% 必然有一个共同的，删去
                                        path_cell{nk}{end+1}=fliplr([path_t(1:end-1);path_cell{qq}{kk}]);% 必然有一个共同的，删去
                                        % 添加邻居
                                        neighbor_cell{nk}(end+1)=pp;
                                        neighbor_cell{pp}(end+1)=nk;

                                        path_cell{nk}(neighbor_cell{nk}==qq)=[];
                                        neighbor_cell{nk}(neighbor_cell{nk}==qq)=[];
                                    end
                                end
                            end
                            ks{pp}=[ks{pp};ks{qq}];
                            pore_c_id{pp}=unique([pore_c_id{pp} pore_c_id{qq}]);% 添加新的相邻的圆柱
                        end
                    end
                end
            end
        end
    end
    %%
    save(['..\basic\pack\' fileList(nn).name '\pore_1.mat'],'ks','throat_cell','path_cell','is_del','neighbor_cell','pore_c_id');
end
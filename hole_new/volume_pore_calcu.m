% 读取圆盘位置
clear
d=30.75*2/0.8;
h=4.81*2/0.8;% 由于测量不准，需要使圆柱厚一点
r=d/2;
cylinder_size=[r,h];%
fileList=dir('..\basic\pack\pack*');
% load('p.mat');p=p';p_length=length(p);
for nn=1:30
    %% 加载数据
    disp(nn)
    load(['..\basic\pack\' fileList(nn).name '\basic.mat'])
    load(['..\basic\pack\' fileList(nn).name '\rc_is_in.mat'])
    load(['..\basic\pack\' fileList(nn).name '\voro.mat'],'v')
    load(['..\basic\pack\' fileList(nn).name '\pore_1.mat'],'ks','is_del')
    %%
    pore_tet=cell(1,length(ks));
    volume=cell(1,length(ks));    
    parfor ii=1:length(ks)
        pore_tet{ii}=zeros(0,4);
        volume{ii}=0;
        if is_del(ii)
            continue
        end
        pore_tri=ks{ii};
        for jj=1:length(pore_tri)
            tet_neighbor=find(sum(ismember(pore_tri(jj+1:end,:),pore_tri(jj,:)),2)==2);
            tet_neighbor=tet_neighbor+jj;
            for kk=1:length(tet_neighbor)
                tet_tmp=sort(unique([pore_tri(tet_neighbor(kk),:),pore_tri(jj,:)]));
                if ~ismember(tet_tmp,pore_tet{ii},'rows')
                    pore_tet{ii}=[pore_tet{ii};tet_tmp];
                end
            end
        end
        for jj=1:length(pore_tet{ii})
             volume{ii}=volume{ii}+det([v(pore_tet{ii}(jj,:),:) ones(4,1)])/6;% 四面体体积公式
        end
        % 先计算四面体体积之和再减去相交的体积
    end
    %%

    %%

end
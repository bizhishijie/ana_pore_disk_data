idx_keep=setxor(1:size(Rc,2),idx_remove);
Rc=[Rc(:,idx_keep) Rc_add];
Ori=[Ori(:,idx_keep) Ori_add];

[~,idx_eff2]=ismember(idx_eff,idx_keep);
idx_eff2=idx_eff2(idx_eff2~=0);
idx_eff=[idx_eff2 length(idx_keep)+1:1:length(idx_keep)+size(Rc_add,2)];

% clearvars -except Rc Ori idx_eff  pack_num load_address D H
function cell_out = cell_combine_ (cell_in)
%% combine index in cell_in that are in different cells


idx_remain=1:length(cell_in);
cell_out=cell(1,0);
while ~isempty(idx_remain)
    idx_start=[];
    idx_start_new=idx_remain(1);
    while length(idx_start_new)~=length(idx_start)
        idx_start=idx_start_new;
        idx_start_add=cell2mat(cell_in(idx_start));
        idx_start_new=unique([idx_start idx_start_add]);
    end
    cell_out{length(cell_out)+1}=idx_start_new;
    idx_remain=setxor(idx_remain,idx_start_new);
end



end
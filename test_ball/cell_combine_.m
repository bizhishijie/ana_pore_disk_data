function cell_out = cell_combine_ (cell_in)
a=1:length(cell_in);
a=repelem(a,cellfun('length',cell_in));
b=cell2mat(cell_in);
g=graph(a,b);
cell_out=conncomp(g,'OutputForm','cell');
end
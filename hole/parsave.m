function parsave(parsave_file_name, varargin)
var_names={};
for i=1:length(varargin)
    eval(sprintf('%s=varargin{%d};',inputname(i+1),i));
    var_names{i}=inputname(i+1);
end
save(parsave_file_name,var_names{:});
end

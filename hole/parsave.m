function parsave(fname,data)
var_name=genvarname(inputname(2));
eval([var_name '=data']);

save(fname,var_name);
end
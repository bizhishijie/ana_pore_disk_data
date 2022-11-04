function is_in=is_in_box(point,box)
min_box=box(:,1);
max_box=box(:,2);

flag1=any(point<min_box');
flag2=any(point>max_box');
is_in=~(flag1&flag2);
end
function is_in=is_point_in_tetra(p,p1,p2,p3,p4)
a=p2-p1;
b=p3-p1;
c=p4-p1;
p=p-p1;

aa=dot(a,a);
ab=dot(a,b);
bb=dot(b,b);
ac=dot(a,c);
bc=dot(b,c);
cc=dot(c,c);

pa=dot(p,a);
pb=dot(p,b);
pc=dot(p,c);

f=[aa ab ac;ab bb bc;ac bc cc]\[pa;pb;pc];
is_in=all(f>=0)&all(f<=1)&sum(f)>=0&sum(f)<=1;
end
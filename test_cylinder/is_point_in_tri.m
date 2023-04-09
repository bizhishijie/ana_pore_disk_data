function is_in=is_point_in_tri(p,p1,p2,p3)
a=p2-p1;
b=p3-p1;
p=p-p1;
ab=dot(a,b);
aa=dot(a,a);
bb=dot(b,b);
pa=dot(p,a);
pb=dot(p,b);
f=[aa ab;ab bb]\[pa;pb];
is_in=f(1)>0&&f(1)<1&&f(2)>0&&f(2)<1&&f(1)+f(2)>0&&f(1)+f(2)<1;
end
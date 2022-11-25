d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
fd=@(p) dcylinder(p,r,h/2,18);
[p,t]=distmeshsurface(fd,@huniform,1.5,1.1*[-r,-r,-h/2;r,r,h/2]);
p=p';t=t';
save('..\hole\tri.mat','p','t')
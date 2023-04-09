function d=dcylinder(p,r,h_2,n)
% d=30.75*2;
% h=4.81*2;
% r=d/2;
% fd=@(p) dcylinder(p,r,h/2,16);
% [p,t]=distmeshsurface(fd,@huniform,0.05,1.1*[-r,-r,-r;r,h/2,h/2]);
% save('..\hole\tri.mat')
d=sqrt(p(:,1).^2+p(:,2).^2).^n/r^n+p(:,3).^n/h_2^n-1;

end
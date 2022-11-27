d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
h_num=3;
hl=linspace(-h/2+h/h_num/2,h/2-h/h_num/2,h_num);
fd=@(p) sqrt(sum(p.^2,2))-r;
fstats=@(p,t) fprintf('%d nodes, %d elements, min quality %.2f\n', ...
                      size(p,1),size(t,1),min(simpqual(p,t)));
[p0,t]=distmesh2d(fd,@huniform,6,[-r -r;r r],[]);
fstats(p0,t);
dis=fd(p0);
p=zeros(length(p0)*2+ sum(dis>-1)*h_num,3);
p(1:length(p0),1:2)=p0;p(1:length(p0),3)=h/2;
p(length(p0)+1:2*length(p0),1:2)=p0;p(length(p0)+1:2*length(p0),3)=-h/2;
p(2*length(p0)+1:end,:)=[repmat(p0(dis>-1,:),h_num,1) sort(repmat(hl,1,sum(dis>-1)))'];
save('..\hole\p.mat','p')
plot3(p(:,1),p(:,2),p(:,3),'o')
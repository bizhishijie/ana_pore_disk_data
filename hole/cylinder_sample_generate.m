d=30.75*2/0.8;
h=4.81*2/0.8;
r=d/2;
r_num=15;
th_num=15;
h_num=5;
ra=linspace(r/(r_num-1),r,r_num);
th=linspace(2*pi/th_num,2*pi,th_num)';
ha=linspace(-h/2+h/2/h_num,h/2-h/2/h_num,h_num);
x=cos(th)*ra;x=x(:);
y=sin(th)*ra;y=y(:);
p=[[x y h/2*ones(length(x),1)]',[x y -h/2*ones(length(x),1)]'];
x=r*cos(th);x=repmat(x,h_num,1);
y=r*sin(th);y=repmat(y,h_num,1);
ha=sort(repmat(ha',th_num,1));
p=[p [x y ha]'];

[~, I, ~] = unique(p,'first','rows');
I = sort(I);
p = p(I,:);

plot3(p(1,:),p(2,:),p(3,:),'o')

save('tri.mat','p')
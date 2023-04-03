close all

alltestcase = -1:8;
ncases = length(alltestcase);
profile on
for i=ncases:-1:1
    testcase = alltestcase(i);
    switch testcase
        case -1
            A0 = [-2 -1 0;
                4 3 0].';
            r0 = 2.2;
            A1 = [10 -1 1;
                16 3 1].';
            r1 = 0.3;
        case 0
            A0 = [-2 -1 0;
                4 3 0].';
            r0 = 0.2;
            A1 = [-2 -1 1;
                4 3 1].';
            r1 = 0.3;
        case 1
            A0 = [-2 -1 0;
                4 3 0].';
            r0 = 0.2;
            A1 = [-1 2 1;
                5 -1 2].';
            r1 = 0.3;
        case 2
            A0 = [-2 -1 0;
                4 3 0].';
            r0 = 0.5;
            A1 = [3 0 1;
                5 -1 0].';
            r1 = 0.3;
        case 3
            A0 = [-2 -1 0;
                1 1 0].';
            r0 = 2;
            A1 = [3 2 1;
                5 1 1].';
            r1 = 1;
        case 4
            A0 = [-2 -1 0;
                1 1 0].';
            r0 = 1;
            A1 = [3 0 1;
                5 -1 0].';
            r1 = 0.5;
        case 5
            A0 = [-2 0 0;
                5 0 0].';
            r0 = 1;
            A1 = [1 3 1;
                1 7 1].';
            r1 = 1.2;
        case 6
            A0 = [-2 0 0;
                5 0 0].';
            r0 = 0.5;
            A1 = [1 2 0;
                  4 1 0].';
            r1 = 0.5;
        case 7
            A0 = [-3 0 0;
                  0 0 0].';
            r0 = 3;
            A1 = [0.5 0 -0.5;
                  3.5 3 0].';
            r1 = 2;
         case 8
            A0 = [-3 0 0;
                  0 0 0].';
            r0 = 0.2;
            A1 = [0.5 0 0.2;
                  -3 -3 0].';
            r1 = 0.2;             
    end
    
    [Distance,P0,P1] = dcylinders(A0, r0, A1, r1);
    s = struct('A0',A0,...
               'r0',r0,...
               'A1',A1,...
               'r1',r1,...
               'Distance',Distance,...
               'P0',P0,...
               'P1',P1,...
               'dummy',[]...
               );
    stab(i) = s;
    Distance
end

profile off
profile viewer

for i=1:ncases
    s = stab(i);
    fig = figure();
    ax = axes('Parent',fig);
    hold(ax,'on');
    PlotCyl(ax, s.A0, s.r0);
    PlotCyl(ax, s.A1, s.r1);
    axis(ax,'equal');
    P = [s.P0 s.P1].';
    plot3(ax,P(:,1),P(:,2),P(:,3),'-r.');
end

%%
function PlotCyl(ax, A, r)

ntt = 16;
theta = linspace(0,2*pi,ntt+1);
circle = [cos(theta); sin(theta)];
cylaxes = A;
D = diff(A,1,2);
Q = null(D.');
c = (r*Q)*circle;
c1 = c + cylaxes(:,1);
c2 = c + cylaxes(:,2);
x = [c1(1,1:end-1); c1(1,2:end); c2(1,2:end); c2(1,1:end-1)];
y = [c1(2,1:end-1); c1(2,2:end); c2(2,2:end); c2(2,1:end-1)];
z = [c1(3,1:end-1); c1(3,2:end); c2(3,2:end); c2(3,1:end-1)];
c = permute(cat(3,c1,c2),[2 3 1]);
%plot3(ax, x, y, z, 'Color', 0.7+[0 0 0]);
%plot3(ax, c(:,:,1), c(:,:,2), c(:,:,3), 'Color', 0.7+[0 0 0])
fill3(ax, x, y, z, 0.7+[0 0 0],'FaceAlpha',0.5);
fill3(ax, c(:,:,1), c(:,:,2), c(:,:,3), 0.7+[0 0 0],'FaceAlpha',0.5)

end

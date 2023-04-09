function [Distance,P0,P1] = dsphereplane(center, rs, p, pnormflag) %#ok
% [Distance,P0,P1] = dsphereplane(Ac, rc, p, true)
% Distance between sphere and a plane in 3D
% INPUTS:
%   center: 3x1 array, sphere center
%   rs: scalar, radius of cylinder 1
%   p:  (1 x 4) equation of the plane, meaning d = p*[x;y;z;1] give a
%       the true distance of (x,y,z) to the plane. In this case one can show
%       norm(p1:3) == 1 
% If p is not normalized, meaning  norm(p1:3) ~= 1, you must call
% dsphereplane(Ac, rc, p) so that dsphereplane will normalize p
% internally
% OUTPUTS:
%   Distance: scalar, if positive, the distance between sphere and plane
%             meaning Distance = min{ |x-y| : x in sphere, y in plane }
%             where |.| is euclidian distance
%   If Distance < 0 the sphere cross the plane
%   P0, P1 closest points, P0 in sphere, P1 in plane such that
%          Distance = |P0-P1|.
% Required: MATLAB 2016b (auto expansion)
%
% Test script
%
% close all
% clear
% ncases = 1;
% for i=ncases:-1:1
%     center = randn(3,1);
%     r = rand();
%     p = rand(1,4);
%     %load toto.mat
%     [Distance,P0,P1] = dsphereplane(center, r, p);
%     s = struct('center',center,...
%                'r',r,...
%                'p', p,...
%                'Distance',Distance,...
%                'P0',P0,...
%                'P1',P1,...
%                'dummy',[]...
%                );
%     stab(i) = s;
%     Distance;
%     if (norm(P0-P1)-abs(Distance))/abs(Distance) > 1e-10
%         keyboard
%     end
%     d1 = p/norm(p(1:3))*[center; 1];
%     if any(d1 < Distance)
%         keyboard
%     end
% end
% 
% 
% [X,Y,Z] = sphere();
% for i=1:ncases
%     s = stab(i);
%     fig = figure();
%     ax = axes('Parent',fig);
%     hold(ax,'on');
%     x = X*s.r+s.center(1);
%     y = Y*s.r+s.center(2);
%     z = Z*s.r+s.center(3);
%     surf(ax,x,y,z);
%     axis(ax,'equal');
%     P = [s.P0 s.P1].';
%     plot3(ax,P(:,1),P(:,2),P(:,3),'-r.');
%     Q = null(s.p(1:3));
%     theta = linspace(0,2*pi);
%     c = [cos(theta); sin(theta)];
%     C = s.P1 + Q*c;
%     patch(ax,C(1,:),C(2,:),C(3,:),'b');
% end
% %%
% function PlotCyl(ax, A, r)
% ntt = 16;
% theta = linspace(0,2*pi,ntt+1);
% circle = [cos(theta); sin(theta)];
% cylaxes = A;
% D = diff(A,1,2);
% Q = null(D.');
% c = (r*Q)*circle;
% c1 = c + cylaxes(:,1);
% c2 = c + cylaxes(:,2);
% x = [c1(1,1:end-1); c1(1,2:end); c2(1,2:end); c2(1,1:end-1)];
% y = [c1(2,1:end-1); c1(2,2:end); c2(2,2:end); c2(2,1:end-1)];
% z = [c1(3,1:end-1); c1(3,2:end); c2(3,2:end); c2(3,1:end-1)];
% c = permute(cat(3,c1,c2),[2 3 1]);
% fill3(ax, x, y, z, 0.7+[0 0 0],'FaceAlpha',0.5);
% fill3(ax, c(:,:,1), c(:,:,2), c(:,:,3), 0.7+[0 0 0],'FaceAlpha',0.5)
% end


N = p(1:3).';
if nargin < 4
    nN = sqrt(N'*N);
    p = p / nN;
    N = N / nN;
end

dc = p*[center;1];
Distance = dc - rs;

P0 = center-N*rs;
P1 = P0 - N*Distance;

end

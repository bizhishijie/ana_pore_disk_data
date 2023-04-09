function [Distance,P0,P1] = dsphere(C0, r0, C1, r1)
% [Distance,P0,P1] = dsphere(C0, r0, C1, r1);
% Distance between 2 spheres in 3D
% INPUTS:
%   C0, C1: 3x1 arrays, centers of the spheres
%   r0, r1: scalar radius of the sphere
% OUTPUTS:
%   Distance: scalar, if positive, the distance between spheres
%             meaning Distance = min{ |x-y| : x in sphere1, y in sphere2 }
%             where |.| is euclidian distance
%   If Distance < 0 the objects intersect
%   P0, P1 closest points, P0 in sphere1, P1 in sphere2 such that
%          Distance = |P0-P1|.
% Required: MATLAB 2016b (auto expansion)
%
% Author: Bruno Luong <brunoluong@yahoo.fr>

D = C1-C0;
Distance = sqrt(D'*D) - (r0+r1);

dP = D/norm(D);
P0 = C0 + dP*r0;
P1 = C1 - dP*r1;

end

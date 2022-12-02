function [Distance,P0,P1] = distance_cylinder_sphere(Ac, rc, C, rs)
% [Distance,P0,P1] = dcylindersphere(Ac, rc, A1, rs)
% Distance between cylinder and a sphere in 3D
% INPUTS:
%   Ac: 3x2 array, each column is the end of axis of cylinder 0
%   rc: scalar, radius of cylinder 0
%   center: 3x1 array, sphere center
%   rs: scalar, radius of cylinder 1
% OUTPUTS:
%   Distance: scalar, if positive, the distance between cylinder and sphere
%             meaning Distance = min{ |x-y| : x in cylinder, y in sphere }
%             where |.| is euclidian distance
%   If Distance < 0 the objects intersect
%   P0, P1 closest points, P0 in cylinder, P1 in sphere such that
%          Distance = |P0-P1|.
% Required: MATLAB 2016b (auto expansion)
%
% Ac=[0 0; 0 0; -3 1]; rc=0.5;
% C=[2; 0; 0]; rs=1;
% [Distance,P0,P1] = dcylindersphere(Ac, rc, C, rs) % return 0.5, [0.5;0;0]; [1;0;0]
%
% Ac=[0 0; 0 0; -1 1]; rc=1;
% C=[4; 0; 5]; rs=1;
% [Distance,P0,P1] = dcylindersphere(Ac, rc, C, rs) % return 4
%
% Author: Bruno Luong <brunoluong@yahoo.fr>

Dc = diff(Ac,1,2);
t = Dc'*(C-Ac(:,1));
Dc2 = Dc'*Dc;
if (t >=0) && (t <= Dc2)
    CA1 = C-Ac(:,1);
    t2 = t*t;
    t2 = t2/Dc2;
    D2 =  CA1'*CA1 - t2;
    Distance = sqrt(max(D2,0)) - (rc+rs);
    P0 = Ac(:,1) + sqrt(t2/Dc2)*Dc;
    dP = C-P0;
    P0 = P0+dP*(rc/norm(dP));
else
    Q = null(Dc');
    if t < 0
        P0 = Ac(:,1);
        Ac(:,1);
    else
        P0 = Ac(:,2);
        t = t-Dc2;
    end
    CA = (C-P0);
    t2 = t*t/Dc2;
    Cp = Q'*CA;
    Cp2 = Cp'*Cp;
    if Cp2 < rc*rc
        Pb = Cp;
        dCcap2 = 0;
    else
        Pb = Cp*(rc/sqrt(Cp2));
        dCcap2 = (sqrt(Cp2)-rc)^2;
    end
    P0 = P0 + Q*Pb;
    Distance = sqrt(dCcap2 + t2) - rs;
end

dP = P0-C;
P1 = C + dP*(rs/norm(dP));

end

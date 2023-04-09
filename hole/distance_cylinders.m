function [Distance,P0,P1] = distance_cylinders(A0, r0, A1, r1)
% [Distance,P0,P1] = dcylinders(A0, r0, A1, r1)
% Distance between 2 cylinders in 3D
% INPUTS:
%   A0: 3x2 array, each column is the end of axis of cylinder 0
%   r0: scalar, radius of cylinder 0
%   A1: 3x2 array, each column is the end of axis of cylinder 1
%   r1: scalar, radius of cylinder 1
% OUTPUTS:
%   Distance: scalar, if positive, the distance between 2 cylinders
%             meaning Distance = min{ |x-y| : x in cylinder0, y in cylinder1 }
%             where |.| is euclidian distance
%   If Distance < 0 the cylinders intersect
%   P0, P1 closest points, P0 in cylinder0, P1 in cylinder1 such that
%          Distance = |P0-P1|.
% Required: MATLAB 2016b (auto expansion)
% Author: Bruno Luong <brunoluong@yahoo.fr>

% Compile mexfile CONV1.mex if it doesn't exist
compileconv1();

% cylinder0 shell x cylinder1 shell
SA = SubspaceAngles(A0, r0, A1, r1, [], 1);
[Distance,P0,P1,facekept] = Case1(A0, r0, A1, r1, SA);
if ~isnan(Distance)
    return
end

%% Face on cylinder shell
SA = SubspaceAngles(A0, r0, A1, r1, SA, 2);
Distance = nan(4,1);
P0 = nan(3,4);
P1 = nan(3,4);
for k=find(facekept)
   [Distance(k),P0(:,k),P1(:,k)] = Case2(A0, r0, A1, r1, SA, k);
   if Distance(k) < 0
       break
   end
end
[Distance,imin] = min(Distance);
if ~isnan(Distance)
    P0 = P0(:,imin);
    P1 = P1(:,imin);
    return
end

%% Faces on faces
SA = SubspaceAngles(A0, r0, A1, r1, SA, 3);
Distance = nan(4,1);
P0 = nan(3,4);
P1 = nan(3,4);

C0 = A0(:,[1 1 2 2]);
R0 = r0 + zeros(1,4);
s0 = [-1 -1 +1 +1];

C1 = A1(:,[1 2 1 2]);
R1 = r1 + zeros(1,4);
s1 = [-1 +1 -1 +1];

dC = sqrt(sum((C0-C1).^2,1));
[dCs, ktab] = sort(dC);
dCmin = dCs(1);
ilast = find(dCs<=dCmin+(r0+r1),1,'last');
ktab = ktab(1:ilast);
for k = ktab
   [Distance(k),P0(:,k),P1(:,k)] = Case3(C0(:,k), R0(k), s0(k), ...
                                         C1(:,k), R1(k), s1(k), SA);
   if Distance(k) < 0
       break
   end
end

[Distance,imin] = min(Distance);
P0 = P0(:,imin);
P1 = P1(:,imin);

end

%%
function crc = GetDiscreteCircle()
% Need in case numerical issue
persistent CRC
if isempty(CRC)
    ntt = 144; % increase ntt for more accuracy
    theta = linspace(0,2*pi,ntt+1);
    theta(end) = [];
    ctt = cos(theta);
    stt = sin(theta);
    CRC = [ctt; stt];
end
crc = CRC;
end

%%
function SA = SubspaceAngles(A0, r0, A1, r1, SA, casenum)
% Initialization few quantities we need over and over again

if casenum <= 1
    D0 = diff(A0,1,2);
    D1 = diff(A1,1,2);
    lax0 = norm(D0);
    lax1 = norm(D1);
    N0 = D0 / lax0;
    N1 = D1 / lax1;
    N0xN1 = cross(N0,N1);    
       
    SA1 = struct(...
            'D0',   D0,...
            'D1',   D1,...
            'lax0', lax0,...
            'lax1', lax1,...
            'N0',   N0,...
            'N1',   N1,...
            'N0xN1',N0xN1);
        
    SA = struct('Case1',SA1);    
        
elseif casenum <= 2   
    % Need for case2
    [Q0,~] = qr(SA.Case1.N0);
    Q0(:,1) = []; % orth basis of cylinder cross section
    [Q1,~] = qr(SA.Case1.N1);
    Q1(:,1) = []; % orth basis of circle plane
    M = Q0.'*Q1;
    [U,S,V] = svd(M);
    s = reshape(diag(S),[1 2]);    
    SA.Case2 = struct(...
            'Q0',   Q0,...
            'Q1',   Q1,...
            'M',    M,...
            'U',    U,...
            'V',    V,...
            's',    s);
        
else %  casenum <= 3
    % Need for case3
    UV1 = SA.Case2.Q1;
    if det([UV1,SA.Case1.N1]) < 0
        UV1(:,2) = -UV1(:,2);
    end
    UV1 = r1*UV1;
    N0xUV1 = cross(SA.Case1.N0,UV1);    
    ParallelPlan = norm(SA.Case1.N0xN1) < 1e-9;
    
    if ~ParallelPlan
        %% Determine if the centers are face2face
        C1 = A1(:,1);
        C0 = A0(:,1);
        D = C1-C0;
        P0 = SA.Case2.Q0 * SA.Case2.U(:,2);
        P1 = SA.Case2.Q1 * SA.Case2.V(:,2);
        P01 = [P0,-P1];
        x = P01 \ D;
        dV = D - P01*x;
        dv2 = dV.'*dV;
        Face2Face = (dv2 < 1e-6*min([r0 r1])^2);
        P0r0 = r0*P0;
        P1r1 = r1*P1;
    else
        Face2Face = false;
        P01 = [];
        P0r0 = [];
        P1r1 = [];
    end
    SA.Case3 = struct(...
            'UV1',  UV1, ...
            'N0xUV1',N0xUV1,...
            'ParallelPlan',ParallelPlan,...
            'Face2Face',Face2Face,...
            'P01',P01,...
            'P0r0',P0r0,...
            'P1r1',P1r1);
end

end % SubspaceAngles

%%
function [Distance,P0,P1,facekept] = Case1(A0, r0, A1, r1, SA)
% Case 1, distance between 2 shells
SA = SA.Case1;
D0 = SA.D0;
D1 = SA.D1;
N = SA.N0xN1;
n = norm(N);
if n < 1e-9
    % parallel cylinders
    % t = projection on common axes
    t0 = D0.'*A0;
    t1 = D0.'*A1;
    s0 = sort(t0);
    s1 = sort(t1);
    smin = max([s0(1) s1(1)]);
    smax = min([s0(2) s1(2)]);
    t = 0.5*(smin+smax);
    if smin <= smax
        w1 = (t-t0(1))/(t0(2)-t0(1));
        w0 = 1-w1;
        P0 = A0*[w0;w1];
        w1 = (t-t1(1))/(t1(2)-t1(1));
        w0 = 1-w1;
        P1 = A1*[w0;w1];
        dP = P1-P0;
        d = norm(dP);
        Distance = d-(r0+r1);
        dP = dP/d;
        P0 = P0 + r0*dP;
        P1 = P1 - r1*dP;
    else
        Distance = NaN;
        P0 = nan(3,1);
        P1 = nan(3,1);
    end
    x = [(t-t0(1))/(t0(2)-t0(1)), ...
         (t-t1(1))/(t1(2)-t1(1))];    
    facekept = [x(1)<=0, x(1)>=1, x(2)<=0, x(2)>=1];
else
    N = N/n;
    % axis not parallel, find cross point in the common subspace
    % that contains both axis drections
    dA10 = A1(:,1)-A0(:,1);
    D = [-D0, D1];
    x = -(D\dA10);
    facekept = [x(1)<=0, x(1)>=1, x(2)<=0, x(2)>=1]; 
    if all(~facekept)
        d = N.'*dA10; % abs(d) == norm(M*[x;1]); with M:=[D,dA10]
        % force the normal point from axe0 to axe1
        if (d < 0)
            N = -N;
            d = -d;
        end
        Distance = d - (r0+r1);
        P0 = A0(:,1) + x(1)*D0 + r0*N;
        P1 = A1(:,1) + x(2)*D1 - r1*N;
    else
        Distance = NaN;
        P0 = nan(3,1);
        P1 = nan(3,1);
    end
end

end % Case1

%% Find real roots of a polynomial
function r = rroots(Poly)
r = roots(Poly);
ir = abs(imag(r));
irmax = max(ir);
keep = ir <= 1e-6*irmax;
if any(keep)
    r = real(r(keep));
else
    % well, we know the polynomial must have at least a real toots
    r = real(r);
end
r = reshape(r, 1, []);
end

%%
function [Distance,P0,P1] = Case2(A0, r0, A1, r1, SA, k)
% Circular Face (cap) on cylinder shell

switch k
    case {1,2}
        C = A0(:,k);
        r = r0;
        Ac = A1;
        D = SA.Case1.D1;
        rc = r1;
        Q = SA.Case2.Q1;
        P = SA.Case2.Q0;
        U = SA.Case2.V;
        V = SA.Case2.U;      
    case {3,4}
        C = A1(:,k-2);
        r = r1;
        Ac = A0;
        D = SA.Case1.D0;
        rc = r0;
        Q = SA.Case2.Q0;
        P = SA.Case2.Q1;
        U = SA.Case2.U;
        V = SA.Case2.V;                
end

% %N = N / norm(N);
% Q = null(diff(Ac,1,2).'); % orth basis of cylinder cross section
% P = null(N.'); % orth basis of circle plane
% [U,S,V] = svd(Q.'*P);

x = Q.'*(Ac(:,1)-C);
sr = r * SA.Case2.s;

A = U.*sr;
y = -x.'*A;
dsr2 = sr(1)*sr(1)-sr(2)*sr(2);
if dsr2==0 && ((y*y.') < 1e-10*min([r0 r1])^2)
    % Concentric circles, we peak any point
    T = [1; 0];
else   
    c = [ dsr2, ...
         -y(2),...
          y(1)];
    Poly = c * [ 0 -2  0 2 0;   % P(t) = (2*t) * (1-t^2)
                -1  0  0 0 1;   % P(t) = (1-t^2) * (1+t^2)
                 0  2  0 2 0];  % P(t) = (2*t) * (1+t^2)
    tprj = rroots(Poly);
    
    if isempty(tprj) % issue with numerical inaccuracy
        % brute force
        T = GetDiscreteCircle();
    else
        % Candidate of the projections
        T = [(1-tprj.^2); 2*tprj] ./ (1+tprj.^2);
    end
end

xprj = A*T;
d2 = sum((xprj-x).^2,1);
[d2,imin] = min(d2);
xprj = xprj(:,imin);

P0 = C + r*P*V*T(:,imin); % == C + P*(M\xprj);
dP = Q*(x - xprj);
P1 = P0 + dP;
d0 = sqrt(d2); %==norm(x - xprj);

t = D.'*[Ac, P1];
if (t(3)>=t(1)) && (t(3)<=t(2))
    d = d0 - rc;
    Distance = d;
    if (d > 0)
        P1 = P0 + (d / d0)*dP; 
    end
else
    Distance = NaN;
    P0 = nan(3,1);
    P1 = nan(3,1);
end

end % Case2

%%
function [Distance,P0,P1] = Case3(C0, r0, s0, C1, r1, s1, SA)
% Distance between 2 caps

l0 = SA.Case1.lax0;
l1 = SA.Case1.lax1;
N0 = s0*SA.Case1.N0;
N1 = s1*SA.Case1.N1;

D   = nan(1,3);
P0  = nan(3,3);
P1  = nan(3,3);
[D(1),P0(:,1),P1(:,1)] = dcircledisk(C0, r0, N0, C1, r1, N1, l1, SA, false);
if ~(D(1) < 0) % not the same as ">=0" test if value is NaN
    [D(2),P1(:,2),P0(:,2)] = dcircledisk(C1, r1, N1, C0, r0, N0, l0, SA, true);
    if ~(D(1) < 0)
        [D(3),P1(:,3),P0(:,3)] = dcircles(C0, r0, N0, C1, r1, N1, SA);
    end
end

[Distance,imin] = min(D);
P0 = P0(:,imin);
P1 = P1(:,imin);

end % Case3

%%
function [Distance,P0,P1] = dcircledisk(C0, r0, N0, C1, r1, N1, l, SA, ReverseFlag) %#ok
% Distance between the 
%   1. circumference of face0 (C0,r0,N0) to the surface and
%   2. face1 (C1,r1,N1)
if ReverseFlag
    P = SA.Case2.Q1;
    V = SA.Case2.V;
else
    P = SA.Case2.Q0;
    V = SA.Case2.U;   
end
% Q = null(N1.'); % orth basis of second disk
% P = null(N0.'); % orth basis of first circle
% M = Q.'*P;
% [~, S, V] = svd(M);
% s = diag(S);

% V(:,2) always corresponds to the smallest singular value
P0 = C0 + (r0*[-1 1]).*(P*V(:,2));
alpha = N1.' * (P0-C1);
[~,imin] = min(abs(alpha));
d = alpha(imin);
P0 = P0(:,imin);
P1 = P0 - d * N1;
if alpha(1)*alpha(2) >= 0
    sqr_dP1C1 = sum((P1-C1).^2,1);
    if sqr_dP1C1 > r1*r1 || (d < -l)
        Distance = NaN;
    else
        Distance = d;
    end
else
    % Check if the two cicles are interlaced
    sn = (alpha(1)+alpha(2)) / (alpha(1)-alpha(2)); % == interp1(alpha,[-1 1],0)
    cs = sqrt(1-sn.^2);
    C = [cs, -cs;
         sn,  sn];
    L = r0*P*V*C; % Vector parallel to the intersect to two planes
    D = diff(L,1,2);
    E = C0 - C1 + L(:,1);
    h1 = (2*r0*cs)^2;
    x = -(D.'*E)/h1;
    Cp1 = E + x*D;
    a2 = sum(Cp1.^2); % square of the distance from C1 to the intersect
    dx2 = r1^2-a2;
    if dx2 >= 0
        dx = sqrt(dx2/h1);
        x = x + dx*[-1 1];
        if x(2) < 0 || x(1) > 1
            Distance = NaN;
        else
            P0 = E + C1 + max(0,x(1))*D;
            P1 = E + C1 + min(1,x(2))*D;            
            Distance = d;
        end
    else
        Distance = NaN;
    end
end

end

%%
function [Distance,P0,P1] = dcircles(C0, r0, N0, C1, r1, N1, SA) %#ok
% [Distance,P0,P1] = dcircles(C0, r0, N0, C1, r1, N1)
% Distance between two circles in 3D
% David Eberly, Geometric Tools, Redmond WA 98052

% UV = null(N1.');
% if det([UV N1]) < 0
%     UV(:,2) = -UV(:,2);
% end

D = C1-C0;
if SA.Case3.ParallelPlan
    % Parallel circles
    N = N0;
    d2 = D.'*D;
    if d2 == 0
        Distance = 0;
        P0 = C0;
        P1 = C1;
    else
        h = D.'*N;
        h2 = h*h;
        a2 = d2 - h2;
        sr = r0+r1;
        if a2 >= sr*sr
            Distance = sqrt((sqrt(a2)-sr)^2+h2);
        else
            Distance = abs(h);
        end
        D = D - h*N;
        D = D / norm(D);
        P0 = C0 + D*r0;
        P1 = C1 - D*r1;
    end
    return
end

if SA.Case3.Face2Face
    % centers are face-to-face wrt intersection line of 2 planes
    x = SA.Case3.P01 \ D;
    s0 = x(1);
    s1 = x(2);
    P0 = C0 + sign(s0)*SA.Case3.P0r0;
    P1 = C1 + sign(s1)*SA.Case3.P1r1;
    tol = 1e-15;
    if abs(s0)*(1+tol)<r0 || abs(s1)*(1+tol)<r1
        Distance = 0;
    else
        Distance = norm(P1-P0);
    end
    return
end

UV1     = SA.Case3.UV1;
N0xres  = [SA.Case3.N0xUV1 cross(N0,D)];
r0sqr   = r0*r0;

a01 = (D.' * UV1);
a0 = a01(1); % == r1*(D.'*U1);
a1 = a01(2); % == r1*(D.'*V1);

ps = N0xres.' * N0xres;
a2 = ps(3,3);  % == sum(N0xD.^2);
a3 = ps(3,1);  % == r1*(N0xD.'*N0xU1);
a4 = ps(3,2);  % == r1*(N0xD.'*N0xV1);
a5 = ps(1,1);  % == r1sqr*sum(N0xU1.^2);
a6 = ps(1,2);  % == r1sqr*(N0xU1.'*N0xV1);
a7 = ps(2,2);  % == r1sqr*sum(N0xV1.^2);

p0 = [a5-a7, 2*a3, a2+a7];
p1 = [2*a6, 2*a4];
p4 = [2*a6, a4, -a6];
p5 = [a7-a5, -a3];
mtwoa0a1 = -2*a1*a0;
tmp0 = [-1, 0, 1];
tmp1 = [a1*a1,0,0] + (a0*a0)*tmp0;                              % [a1*a1,0,0] == conv1(p2,p2)
tmp2 = [mtwoa0a1 0];                                            % == p2*(2*p3);
tmp3 = conv1(p4,p4) + conv1(tmp0, conv1(p5,p5));
p6 = conv1(p0,tmp1) + conv1([-p1,p1],tmp2) - r0sqr*tmp3;        % [-p1,p1] == conv1(tmp0,p1)
p7 = [p0*mtwoa0a1,0] + conv1(p1,tmp1) - r0sqr*conv1(p4,2*p5);   % [p0*mtwoa0a1,0] == conv1(p0,tmp2)

if max(abs(p7)) >  max(abs(p6))*1e-10
    phi = conv1(p6,p6) - conv1(tmp0,conv1(p7,p7));
    cs = rroots(phi);
else
    cs = rroots(p6);
    p7 = 0;
end

if isempty(cs) % issue with numerical inaccuracy
    [Distance,P0,P1] = BruteforceC2C(C0, r0, C1, r1, SA);
    return
end
temp = polyval(p7,cs);
sn = -polyval(p6,cs) ./ temp;
isdegen = temp==0 | abs(cs)>1;
c1 = [cs;
      sn];
if any(isdegen)
    csd = min(max(cs(isdegen),-1),1);
    snd = sqrt(1-csd.^2);
    c1(:,isdegen) = [csd; snd];
    % duplicate
    idup = isdegen;
    idup(idup) = (snd~=0);
    c1dup = c1(:,idup);
    c1dup(2,:) = -c1dup(2,:);
    c1 = [c1, c1dup];
end
if any(abs(sum(c1.^2,1)-1) > 1e-3)
    [Distance,P0,P1] = BruteforceC2C(C0, r0, C1, r1, SA);
    return
end

delta = D + UV1 * c1;
N0dDelta = N0.' * delta;
N0dDelta2 = N0dDelta.*N0dDelta;
lenN0xDelta = sqrt(max(sum(delta.^2)-N0dDelta2,0));
diff = lenN0xDelta - r0;
d2 = N0dDelta2 + diff.*diff;
[d2,imin] = min(d2);
Distance = sqrt(d2);
delta = delta(:,imin);
P1 = delta + C0;
delta = delta - N0dDelta(imin).*N0;
P0 = C0 + (r0/sqrt(sum(delta.^2,1))) * delta;

end % dcircles

%%
function [Distance,P0,P1] = BruteforceC2C(C0, r0, C1, r1, SA)
c = GetDiscreteCircle();
UV1 = r1*SA.Case2.Q1;
c1 = C1 + UV1*c;
UV0 = r0*SA.Case2.Q0;
c0 = C0 + UV0*c;
d = reshape(c0.',[],1,3) - reshape(c1.',1,[],3);
d2 = sum(d.^2,3);
[d,imin] = min(d2(:));
Distance = sqrt(d);
[i0,i1] = ind2sub(size(d2),imin);
P0 = c0(:,i0);
P1 = c1(:,i1);
end

%%
function c = cross(a,b)
% Calculate cross product with auto-expansion
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:);
     a(3,:).*b(1,:)-a(1,:).*b(3,:);
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
end

%%
function y = polyval(p,x)
nc = length(p);
% Use Horner's method for general case where X is an array.
y = zeros(size(x));
y(:) = p(1);
for i = 2:nc
    y = x .* y + p(i);
end
end

%%
function R2018a_mexopts = GetR2018a_mexopts()
% R2018a_mexopts = GetR2018a_mexopts()
if ~verLessThan('MATLAB','9.4')
    R2018a_mexopts = {'-R2018a'};
else
    R2018a_mexopts = {};
end
end

%%
function compileconv1()
% compileconv1
if exist('conv1', 'file')~=3
    here = fileparts(mfilename('fullpath'));
    path = cd(here);
    fprintf('Compile conv1');
    R2018a_mexopts = GetR2018a_mexopts();
    mex('-O', R2018a_mexopts{:}, 'conv1.c');
    cd(path);
end
end


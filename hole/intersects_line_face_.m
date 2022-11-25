function rc_inter = intersects_line_face_ (rc_line,rc_face)
% rc_line: coordinates of the line segment
% rc_face: coordinates of face nod

%%
n_face=cross(rc_face(:,2)-rc_face(:,1),rc_face(:,3)-rc_face(:,2));
n_face=n_face/norm(n_face);
e_line=rc_line(:,2)-rc_line(:,1);
e_line=e_line/norm(e_line);
t=-dot(rc_line(:,1)-rc_face(:,1),n_face)/dot(e_line,n_face);
rc_inter=rc_line(:,1)+e_line*t;
if t>norm(rc_line(:,2)-rc_line(:,1))||t<0
    rc_inter=[];
    return
end

e2=rc_face(:,2)-rc_face(:,1);
e3=rc_face(:,3)-rc_face(:,1);

r_1c=rc_inter-rc_face(:,1);
f=[dot(e2,e2) dot(e2,e3);dot(e2,e3) dot(e3,e3)]\[dot(r_1c,e2);dot(r_1c,e3)];

%if sum(f>=0&f<=1)~=2
if ~(f(1)>=0&&f(1)<=1&&f(2)>=0&&f(2)<=1&&sum(f)>=0&&sum(f)<=1)
    rc_inter=[];
end

end


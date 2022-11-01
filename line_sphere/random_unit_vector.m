function vector=random_unit_vector
u=rand;
v=rand;
theta = 2*pi*u;
phi = acos(2*v - 1);
x=sin(theta)*sin(phi);
y=cos(theta)*sin(phi);
z=cos(phi);
vector=[x;y;z];
end
function result = point_in_tetrahedron(point, tetra)
% point: a matrix of size n-by-3, where each row is a point coordinate
% tetra: a matrix of size m-by-12, where each row is a tetrahedron vertex coordinate
% result: a matrix of size n-by-m, where each element is a boolean value
% reshape point to get a matrix of size 1-by-n-by-3
point = reshape(point, 1, [], 3);
% reshape tetra to get four vertices of each tetrahedron
v1 = reshape(tetra(:, 1:3), [], 1, 3);
v2 = reshape(tetra(:, 4:6), [], 1, 3);
v3 = reshape(tetra(:, 7:9), [], 1, 3);
v4 = reshape(tetra(:, 10:12), [], 1, 3);

% calculate the signed volume of each tetrahedron
vol = dot(cross(v2 - v1, v3 - v1, 3), v4 - v1, 3) / 6;

% calculate the signed volume of four sub-tetrahedrons formed by point and three vertices of tetrahedron
vol1 = dot(cross(v2 - point, v3 - point, 3), v4 - point, 3) / 6;
vol2 = dot(cross(v3 - point, v1 - point, 3), v4 - point, 3) / 6;
vol3 = dot(cross(v1 - point, v2 - point, 3), v4 - point, 3) / 6;
vol4 = dot(cross(v1 - point, v3 - point, 3), v2 - point, 3) / 6;

% check if the sign of vol matches the sign of vol1 + vol2 + vol3 + vol4
result = vol == vol1 + vol2 + vol3 + vol4;

end
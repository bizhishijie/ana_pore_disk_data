function [in,on] = inpolygon3d(polygon,point)
% Check if a point is within a planar polygon in threedimensional space. The procedure is to appropriately create projections of the polygon and the point and to apply the Matlab function "inpolygon"
% polygon: Enter as 3-by-N Matrix. There must be 3 rows, one for each coordinate. Number of columns N corresponds to the number of corners of the polygon
% point: Enter as 3-by-1 Vector.
% E-mail: dimitrij.chudinzow@gmail.com
% Check number of corners
number_polygon_corners=size(polygon,2);
if number_polygon_corners<3 % Check if number of corners in polygon is sufficient
    disp('Error: not enough corners in polygon (must be 3 at least)');
else
    % Calculations
    polygon_vector_1=polygon(:,1)-polygon(:,2); % calculating first direction vector of polygon
    polygon_vector_2=polygon(:,1)-polygon(:,3); % calculating second direction vector of polygon
    polygon_normal_vector=cross(polygon_vector_1,polygon_vector_2); % normal vector standing orthogonal on oplygon plane

    vector_polygon_corner_point_Q=polygon(:,1)-point;
    x_axis=[1;0;0];
    y_axis=[0;1;0];
    

    if dot(polygon_normal_vector,vector_polygon_corner_point_Q)~=0 % if point point is not within the plane of the polygon, the dot product of polygon_normal_vector and vector_polygon_corner_point_Q will not be zero
        in=false;
        on=false;

    else % point is within the area of the polygon

        % first we check whether the normal vector of the polygon is parallel to XY-plane or XZ-plane
        crossproduct_polygon_x_axis=cross(polygon_normal_vector,x_axis);
        crossproduct_polygon_y_axis=cross(polygon_normal_vector,y_axis);

        if norm(crossproduct_polygon_y_axis)==0 % polygon is parallel to XZ-plane --> set Y-coordinates to Zero

            polygon_projection=zeros(size(polygon));
            polygon_projection(1,:)=polygon(1,:);
            polygon_projection(3,:)=polygon(3,:);
            polygon_projection(2,:)=0;

            Q_projection=zeros(size(point));
            Q_projection(1)=point(1);
            Q_projection(3)=point(3);
            Q_projection(2)=0;

            [in,on]=inpolygon(Q_projection(1),Q_projection(3),polygon_projection(1,:),polygon_projection(3,:));

        elseif norm(crossproduct_polygon_x_axis)==0 % polygon is parallel to YZ-plane --> set X-coordinates to Zero

            polygon_projection=zeros(size(polygon));
            polygon_projection(2,:)=polygon(2,:);
            polygon_projection(3,:)=polygon(3,:);
            polygon_projection(1,:)=0;

            Q_projection=zeros(size(point));
            Q_projection(2)=point(2);
            Q_projection(3)=point(3);
            Q_projection(1)=0;

            [in,on]=inpolygon(Q_projection(2),Q_projection(3),polygon_projection(2,:),polygon_projection(3,:));

        else % Polygon is parallel to XY-plane or arbitralily tilted --> set Z-coordinates to Zero

            polygon_projection=zeros(size(polygon));
            polygon_projection(1,:)=polygon(1,:);
            polygon_projection(2,:)=polygon(2,:);
            polygon_projection(3,:)=0;

            Q_projection=zeros(size(point));
            Q_projection(1)=point(1);
            Q_projection(2)=point(2);
            Q_projection(3)=0;

            [in,on]=inpolygon(Q_projection(1),Q_projection(2),polygon_projection(1,:),polygon_projection(2,:));

        end
    end
end
end
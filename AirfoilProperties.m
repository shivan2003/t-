function [A, Sy, Sz, Iy, Iz, Iyz] = AirfoilProperties(airfoil)
%SECTIONPROPERTIES Takes an array of airfoil points and returns airfoil,
%centroid, and moments of inertia
n = length(airfoil) + 1;
y = airfoil(1, :);
z = airfoil(2, :);

y(end+1) = y(1);
z(end+1) = z(1);

A = 0;
Sy = 0;
Sz = 0;

Iz = 0;
Iy = 0;
Iyz = 0;

for i = 1:n-1
    % Formulae from https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    A = A + y(i)*z(i+1) - y(i+1)*z(i);
    Sy = Sy + (y(i) + y(i+1)) * ((y(i)*z(i+1) - y(i+1)*z(i)));
    Sz = Sz + (z(i) + z(i+1)) * ((y(i)*z(i+1) - y(i+1)*z(i)));

    % Formulae from https://en.wikipedia.org/wiki/Second_moment_of_area
    Iz = Iz + (y(i)*z(i+1) - y(i+1)*z(i)) * (y(i)^2 + y(i)*y(i+1) + y(i+1)^2);
    Iy = Iy + (y(i)*z(i+1) - y(i+1)*z(i)) * (z(i)^2 + z(i)*z(i+1) + z(i+1)^2);
    Iyz = Iyz + (y(i)*z(i+1) - y(i+1)*z(i)) * (y(i)*z(i+1) + 2*y(i)*z(i) + 2*y(i+1)*z(i+1) + y(i+1)*z(i)); % Moments about y and z axes
end

A = 1/2 * A;
Sy = 1/(6*A) * Sy;
Sz = 1/(6*A) * Sz;

Iz = 1/12 * Iz;
Iy = 1/12 * Iy;
Iyz = 1/24 * Iyz;

% Use parallel axis theorem to move axes to elastic centre

Iz = Iz - A*Sy^2;
Iy = Iy - A*Sz^2;
Iyz = Iyz - A*Sz*Sy;
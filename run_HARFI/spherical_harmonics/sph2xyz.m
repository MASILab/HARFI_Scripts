function [x, y, z] = sph2xyz(theta, phi, r)

% Transform spherical to Cartesian coordinates.
% theta is the angle from z axis, phi is the counterclockwise angle 
% in the xy plane measured from the positive x axis. theta and
% phi must be in radians

[x, y, z] = sph2cart(phi, pi/2 - theta, r);
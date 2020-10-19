function ang = angleBetw2vect2(x1,y1,z1,x2,y2,z2,halfPiFlag)

% calculate the angle between the two vectors difined by (theta1, phi1) and
% (theta2, phi2), respectively. output ranges [0,pi] when halfPiFlag = 0
% (default).
% if theta1, phi1 are sclar values, theta2, phi2 can be matrix of the same
% size, then output is the a matrix of the same size. if all the inputs are
% vectors of same length, then output will be vector of same length,
% ang(n)=angle btwn [theta1(n), phi1(n)] and [theta2(n), phi2(n)]
%
% by XH
% last modified 08/06/08


if ~exist('halfPiFlag', 'var') || isempty(halfPiFlag)
    halfPiFlag = 0;
end

%[x1, y1, z1] = sph2xyz(theta1, phi1, 1);

%[x2, y2, z2] = sph2xyz(theta2, phi2, 1);

chord = sqrt((x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2);
ang = 2 * real(asin(chord/2)); % asin(x) gives complex results if |x|>1

if (halfPiFlag == 1) 
    ang(ang > pi/2) = pi - ang(ang>pi/2);
end
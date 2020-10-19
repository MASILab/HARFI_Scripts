%% Search the peaks of a spherical function represented by SH coefficients
%% using gradient-based numerical optimization technique (Quasi-Newton)
%% 
%%    
%% Example:
%%   [angles] = find_peaks(coeff,init_pos);
%%  where 'init_pos' is an nx2 array giving n initial directions in
%%  pair of [polar, azimuth], e.g. [80 10; 80 110; 80 150] 

function [angles_out, peaks] = find_peaks_odd(coeff, init_pos, degree)
%%=============================================================
%% Project:   Spherical Harmonics
%% Module:    $RCSfile: find_peaks.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2007/12/27 06:23:35 $
%% Version:   $Revision: 1.3 $
%%=============================================================


n_try = size(init_pos,1);
[n,k] = size(coeff);
angles = zeros(n, n_try*2);
peaks = zeros(n,n_try);
%assume dl=2;
%degree = (sqrt((8*k+1))-3)/2;

%options = optimset('GradObj','on','Display', 'off','LargeScale','off', 'MaxIter', 100,  'TolFun',1e-06, 'TolX',1e-08,'Display','iter');
options = optimset('GradObj','on','Display', 'off','LargeScale','off', 'MaxIter', 100,  'TolFun',1e-06, 'TolX',1e-08);

%tic
for i=1:n
     %disp(sprintf(' i = %d',i));
     angles_i = zeros(1, n_try*2);
     for ii=1:n_try
			 [peak_angle] = fminunc(@(x)real_spherical_harmonics_neg(x,coeff(i,:), degree, 1), init_pos(ii,1:2)*pi/180,  options);        
			 angles_i(2*ii-1:2*ii) = peak_angle;
             peaks(i,ii) = real_spherical_harmonics([peak_angle(1) peak_angle(2)],coeff(i,:), degree, 1);
     end
     %angles_i*180/pi
     angles(i,:) = angles_i;
end    
%toc
angles = angles*180/pi;
angles_out = reshape(angles,[2 60])';


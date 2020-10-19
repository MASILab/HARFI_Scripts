function [Y] = compute_Ylm_all(degree,theta,phi,dl,real_or_complex)
% Evaluate truncated spherical harmonic series of degree "degree"
% // at (theta,phi) where theta and phi represent colatitude and longitude, respectively. 
% // theta is polar angle, 
% INPUTS:
%
% degree                - maximum degree of spherical harmonics
% theta and phi         - location to evaluate
% dl                    - {1} for full band; 2 for even order only  
% real_or_complex       - [{'real'} or 'complex'] basis functions.  

if (nargin<3) 
    error('Usage: [basis] = construct_SH_basis(degree(even), points[nx3], [{real}_or_complex]);');  
end
if (nargin<4) 
    dl = 1;  
end
if (nargin<5) 
    real_or_complex = 'real';  
end

if dl==1
    k =(degree + 1)*(degree + 1);  %% k is the # of coefficients
else
  if (dl==2)
  	k =(degree + 2)*(degree + 1)/2;  %% k is the # of coefficients
  else
    error('dl can only be 1 or 2');
  end
end
Y = zeros(k,1);

for l = 0:dl:degree   %%% even order only due to the antipodal symmetry!
    % calculate the spherical harmonics
    Pm = legendre(l,cos(theta')); % legendre part
    Pm = Pm';
    lconstant = sqrt((2*l + 1)/(4*pi));
    if (dl==2) 
        center = (l+1)*(l+2)/2 - l;
    else
        if (dl==1) 
            center = (l+1)*(l+1) - l;
        end
    end
    Y(center) = lconstant*Pm(1);
    for m=1:l
        precoeff = lconstant * sqrt(factorial(l - m)/factorial(l + m));

        switch lower(real_or_complex)
            case 'real'
                 if mod(m,2) == 1
                     precoeff = -precoeff;
                 end
                Y(center + m) = sqrt(2)*precoeff*Pm(m+1).*cos(m*phi);
                Y(center - m) = sqrt(2)*precoeff*Pm(m+1).*sin(m*phi);
            case 'complex'
                if mod(m,2) == 1
                    Y(center + m) = precoeff*Pm(m+1).*exp(i*m*phi);
                    Y(center - m) = -precoeff*Pm(m+1).*exp(-i*m*phi);
                else
                    Y(center + m) = precoeff*Pm(m+1).*exp(i*m*phi);
                    Y(center - m) = precoeff*Pm(m+1).*exp(-i*m*phi);
                end
            otherwise
                error('The last argument must be either \"real\" (default) or \"complex\".');
        end
    end
end

end



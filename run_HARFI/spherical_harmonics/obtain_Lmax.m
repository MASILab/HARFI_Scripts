function [Lmax Nmin] = obtain_Lmax(DirSignal)
% This function computes the maximum order for the spherical harmonic decomposition (Lmax)
% by taking into account the dimensionality of the signal (number of discrete points)
% The above relationship is based on the real spherical harmonic expansion
% Nmin is the minimum number of points necessary to use Lmax

% Erick Jorge Canales-Rodríguez & Lester Melie-García
%
%%%%%%%
pos = ismember(DirSignal, [0 0 0], 'rows');
ind = find(pos);
DirSignal(ind,:) = [];

N = length(DirSignal);
Lmax = 0;
while (Lmax+1)*(Lmax+2)/2 <= N
  Lmax = Lmax + 2;
end
Lmax = Lmax - 2;
Nmin = sum(1:4:2*Lmax+1);

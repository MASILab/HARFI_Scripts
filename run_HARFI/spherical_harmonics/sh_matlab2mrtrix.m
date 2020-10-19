function sh_mrtrix = sh_matlab2mrtrix(sh_matlab,lmax)
% Converts matlab spherical harmonics to mrtrix format
%
% Only flip some signs. Please read sh_mtrix_vs_matlab.txt
%
% The mrtrix toolbox does not use the (-1)^m factor that Bing's toolbox 
% uses:
% "Note that the spherical harmonics equations used here differ slightly from "
% "those conventionally used, in that the (-1)^m factor has been omitted. This "
% "should be taken into account in all subsequent calculations. \n"

% This results in the ODF appearing flipped in the z-dimension. Just write a for loop and do something like:

% Whatever was made negative by matlab, just undo.
% sh_mrtrix = sh_matlab;
% for l = 0:2:lmax
%     center = (l+1)*(l+2)/2 - l;
%     for m=1:l
%         if mod(m,2) == 1
%             sh_mrtrix(center+m) = -sh_mrtrix(center+m);
%             sh_mrtrix(center-m) = -sh_mrtrix(center-m);
%         end
%     end
% end

% if need to flip along x, flip sign of real (flip positive m)
% if need to flip along y, flip sign of complex (flip negative m)

sh_mrtrix = sh_matlab;
for l = 0:2:lmax
    center = (l+1)*(l+2)/2 - l;
    for m=1:l
%         if mod(m,2) == 1
%             sh_mrtrix(center+m) = -sh_mrtrix(center+m);
%             sh_mrtrix(center-m) = -sh_mrtrix(center-m);
%         end
        %sh_mrtrix(center+m) = -sh_mrtrix(center+m);
        sh_mrtrix(center-m) = -sh_mrtrix(center-m);
    end
end


end





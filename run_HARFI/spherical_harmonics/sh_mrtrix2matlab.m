function sh_matlab = sh_mrtrix2matlab(sh_mrtrix,lmax)
    % Converts mrtrix spherical harmonics to matlab format
    % 
    % Only flip some signs. Please read sh_mtrix_vs_matlab.txt
    % 
    % Note: lmax = -(3/2) + sqrt((3/2)^2 - 4 * (1/2) * (1-length(sh_matlab)))
    % But require lmax anyway since the order should usually always be
    % known
    
    % Whatever was made negative by matlab, just undo.
    sh_matlab = sh_mrtrix;
    for l = 0:2:lmax
        center = (l+1)*(l+2)/2 - l;
        for m=1:l
            if mod(m,2) == 1
                sh_matlab(center+m) = -sh_matlab(center+m);
                sh_matlab(center-m) = -sh_matlab(center-m);
            end
        end
    end
end













function sh_camino = sh_matlab2camino(sh_matlab,lmax)
    % Converts matlab spherical harmonics to camino format
    % 
    % Must reorder AND rescale. Please read sh_matlab_vs_camino.txt
    % 
    % To handle scaling:
    %     camino_scale * camino_coef * Y_unscaled = matlab_scale * matlab_coef * Y_unscaled
    %     => camino_coef = (matlab_scale * matlab_coef)/camino_scale
    %
    % Note: lmax = -(3/2) + sqrt((3/2)^2 - 4 * (1/2) * (1-length(sh_matlab)))
    % But require lmax anyway since the order should usually always be
    % known
    
    % Apply matlab scale
    sh_matlab_scaled = zeros(size(sh_matlab));
    for l = 0:2:lmax
        center = (l+1)*(l+2)/2 - l;
        sh_matlab_scaled(center) = sh_matlab(center);
        for m=1:l
            matlab_scale = sqrt(2);
            if mod(m,2) == 1
                matlab_scale = -matlab_scale;
            end
            sh_matlab_scaled(center+m) = sh_matlab(center+m)*matlab_scale;
            sh_matlab_scaled(center-m) = sh_matlab(center-m)*matlab_scale;
        end
    end
    
    % Reorder from matlab to camino format
    camino_idx = zeros(size(sh_matlab));
    counter = 0;
    for l = 0:2:lmax
        prev_counter = counter;
        for m = 0:l
            counter = counter + 1;
            camino_idx(prev_counter + (2*l+1) - 2*m) = counter;
        end
        for m = 1:l
            counter = counter + 1;
            camino_idx(prev_counter + 2*m) = counter;
        end
    end
    sh_camino_scaled = sh_matlab_scaled(camino_idx);
    
    % Divide by camino scale  
    sh_camino = zeros(size(sh_matlab)); 
    counter = 1;
    for l = 0:2:lmax
        sh_camino(counter) = sh_camino_scaled(counter);
        counter = counter+1;
        for m = 1:l
            sh_camino(counter) = sh_camino_scaled(counter)/2;
            sh_camino(counter+1) = sh_camino_scaled(counter+1)/(-2);
            counter = counter + 2;
        end
    end
end
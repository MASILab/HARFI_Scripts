function sh_matlab = sh_camino2matlab(sh_camino,lmax)
% see sh_matlab2camino for smart sounding stuff
% also run the following lines if you want more info:
% x1 = [119 101 98 40 39 104 116 116 112 58 47 47 119 119 119];
% x2 = [ 46 121 111 117 116 117 98 101 46 99 111 109 47 119 97];
% x3 = [116 99 104 63 118 61 100 81 119 52 119 57 87 103 88];
% x4 = [ 99 81 39 44 39 45 98 114 111 119 115 101 114 39 41];
% eval(char([x1 x2 x3 x4]))

% MULTIPLY by camino scale
sh_camino_scaled = zeros(size(sh_camino));
counter = 1;
for l = 0:2:lmax
    sh_camino_scaled(counter) = sh_camino(counter);
    counter = counter + 1;
    for m = 1:l
        sh_camino_scaled(counter) = sh_camino(counter)*2;
        sh_camino_scaled(counter+1) = sh_camino(counter+1)*(-2);
        counter = counter + 2;
    end
end

% Reorder from matlab to camino format
   matlab_idx = zeros(size(sh_camino));
counter = 0;
for l=0:2:lmax
    prev_counter = counter;
    for m=0:l
        counter = counter+1;
        matlab_idx(counter) = prev_counter + (2*l+1) - 2*m;
    end
    for m=1:l
        counter = counter + 1;
        matlab_idx(counter) = prev_counter + 2*m;
    end
end
sh_matlab_scaled = sh_camino_scaled(matlab_idx);


%% de-apply matlab scale
sh_matlab = zeros(size(sh_camino));
for l=0:2:lmax
    center = (l+1)*(l+2)/2 - l;
    sh_matlab(center) = sh_matlab_scaled(center);
    for m=1:l
        matlab_scale = sqrt(2);
        if mod(m,2) == 1
            matlab_scale = -matlab_scale;
        end
        sh_matlab(center+m) = sh_matlab_scaled(center+m)/matlab_scale;
        sh_matlab(center-m) = sh_matlab_scaled(center-m)/matlab_scale;
    end
end

   
    
   
end



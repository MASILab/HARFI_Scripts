function []= showSH(SH,lmax,spherenum,complex_or_real,even_or_odd)
% SH must be column vector of real_or_complex, 
% even_or_odd: 1 for odd, 2 for even order SH
    

if (nargin<3) 
    spherenum = 60;  
end

% V needs to be Nby1
if any(size(SH)==1)
    if size(SH,2) ~=1
        SH = SH';
    end
else
    error('SH must be column vector')
end


[sphere_x,sphere_y,sphere_z] = sphere(spherenum);
[sh_sphere,~,~] = construct_SH_basis(lmax, [sphere_x(:) sphere_y(:) sphere_z(:)], even_or_odd, complex_or_real);
% 2 designates even, 'real' designates real valued SH

% Multiply basis by sh coefficients to get values of OD on sphere
OD = reshape(sh_sphere * SH,[spherenum+1, spherenum+1]);

% cannot be comple
OD = real(OD);
%OD = abs(OD);

% Remove negative values
OD(OD < 0) = 0;

OD_surf_x = sphere_x.*OD;  % x and y are flipped
OD_surf_y = sphere_y.*OD;
OD_surf_z = sphere_z.*OD;

OD_surf_norm = sqrt(OD_surf_x.^2 + OD_surf_y.^2 + OD_surf_z.^2);
OD_color = cat(3,abs(OD_surf_x./OD_surf_norm),abs(OD_surf_y./OD_surf_norm),abs(OD_surf_z./OD_surf_norm));
OD_color = OD_color(:,:,[1 2 3]); % only here in case we need to flip colors

% Plot as surface
% f = figure();
% parent_axes = axes('parent',f);

surf_h = surf(OD_surf_x, OD_surf_y,OD_surf_z, OD_color);

set(surf_h, 'FaceLighting', 'gouraud', 'FaceColor', 'interp', ...
            'AmbientStrength', 0.8, 'FaceAlpha', 1)
shading(gca,'interp');    
axis equal;

light_h = light('Position', [1, -1, 1], 'Style', 'infinite');
light_h = light('Position', [-1, 1, -1], 'Style', 'infinite');
grid off

xlabel('X')
ylabel('Y')
zlabel('Z')

whitebg(gcf,[0 0 0]);         % background color (in+out plot)
set(gcf,'Color',[0,0,0]);



function [] = plotODsurface(OD_surf_x,OD_surf_y,OD_surf_z,alpha,color)
% plot surface given values along x, y, and z

OD_surf_norm = sqrt(OD_surf_x.^2 + OD_surf_y.^2 + OD_surf_z.^2);

switch nargin
    case 4
        OD_color = cat(3,abs(OD_surf_x./OD_surf_norm),abs(OD_surf_y./OD_surf_norm),abs(OD_surf_z./OD_surf_norm));
        OD_color = OD_color(:,:,[1 2 3]); % only here in case we need to flip colors
    case 5
        OD_color = zeros([size(OD_surf_x,1) size(OD_surf_x,2) 3]);
        OD_color(:,:,1) = color(1);
        OD_color(:,:,2) = color(2);
        OD_color(:,:,3) = color(3);
end

% Plot as surface
% f = figure();
% parent_axes = axes('parent',f);

surf_h = surf(OD_surf_x, OD_surf_y,OD_surf_z, OD_color);

set(surf_h, 'FaceLighting', 'gouraud', 'FaceColor', 'interp', ...
            'AmbientStrength', 0.8, 'FaceAlpha', alpha)
shading(gca,'interp');    
axis equal;
delete(findall(gcf,'Type','light'))
light_h = light('Position', [1, -1, 1], 'Style', 'infinite');
light_h = light('Position', [-1, 1, -1], 'Style', 'infinite');
grid off

xlabel('X')
ylabel('Y')
zlabel('Z')

whitebg(gcf,[0 0 0]);         % background color (in+out plot)
set(gcf,'Color',[0,0,0]);


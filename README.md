# HARFI_Scripts
HARFI processing scripts for MATLAB

calcHARFI.m: runs HARFI on resting state (or task) fMRI data. Returns .mat file with voxel-wise correlation coefficients over spheres

fitHARFI.m: runs HARFI fitting on correlation coefficient data. Returns Spherical Harmonic fits for each voxel, which can be displayed using diffusion MRI [visualization tools](https://github.com/MASILab/dwmri_visualizer)

Jones60.mat: directions over sphere

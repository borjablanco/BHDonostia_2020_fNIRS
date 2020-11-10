function [Jmask, Jmask2] = Jacobian_display_masking_speech (GMSurfaceMesh, J_GM_wav1)
%load('/Users/borjablanco/Documents/MATLAB/my_scripts/RS_reconstruction/J_GM_wav1.mat')
%load('/Users/borjablanco/Documents/MATLAB/head_models/SixMonthHeadModel_updated/Lower_density/GMSurfaceMesh_6months.mat')
clc

% Remove occipital
Jsum = sum(J_GM_wav1);
JsumNorm = abs(Jsum)./max(abs(Jsum(:)));

% Way of displaying the Jacobian
figure
displayIntensityOnMesh_RJC(GMSurfaceMesh, log10(JsumNorm));
caxis([-3 0])

% Masking for group ICA plots
Jmask = JsumNorm>0.01;
Jmask2 = JsumNorm>0.05;
figure
displayIntensityOnMesh_RJC(GMSurfaceMesh, log10(double(Jmask)));
caxis([-1 0])
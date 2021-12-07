% clear; close all; clc

function gPPI_visualization(betas, seed)

% Load required files
addpath(genpath('/Users/borjablanco/Documents/MATLAB/my_scripts/Omniscripts_20190510'))

load('/Users/borjablanco/Documents/MATLAB/my_scripts/borja/greyJet.mat')
load ('/Users/borjablanco/Documents/GitHub/BHDonostia_2020_fNIRS/data_files/GMSurfaceMesh_6months.mat');
load('/Users/borjablanco/Documents/MATLAB/my_scripts/SPEECH_reconstruction/J_GM_wav1.mat')

% Load SD file
load('/Users/borjablanco/Documents/BCBL/Speech_data/Preprocessing/gamma_adjusted/NIRS_SPEECH_5231_preprocessed.mat')

% Load S-D and ch positions
load ('/Users/borjablanco/Documents/MATLAB/my_scripts/SPEECH_reconstruction/ch_pos_GM.mat')
load('/Users/borjablanco/Documents/MATLAB/my_scripts/SPEECH_reconstruction/source_pos.mat')
load('/Users/borjablanco/Documents/MATLAB/my_scripts/SPEECH_reconstruction/detector_pos.mat')

% Load channel position in GM (homer order)
Mpos = ch_pos_GM;
ch = 24;

% Plot GM surface and channels in GM
%figure
%hgm = trisurf(GMSurfaceMesh.face,GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3));hold on
%set(hgm,'edgealpha',0.05,'facecolor','interp','facelighting','phong'),hold on
%alpha(hgm,0.5)
%plot3(ch_pos_GM(:,1),ch_pos_GM(:,2),ch_pos_GM(:,3),'ko','markerfacecolor','k','markersize',6);
%xlabel('x [mm]'),ylabel('y [mm]'),zlabel('z [mm]'),daspect([1 1 1]);grid on
% Plot channel number
%for i = 1:ch
%    shift = 2*[-2*(Mpos(i,1)< mean(Mpos(:,1))) 1 0];
%    text('position',Mpos(i,:)+shift,'string',num2str(i),'fontweight','Bold','fontsize',14,'color','r');
%end

% Smooth GM surface for interpolation
smooth_value = 20;
GMSurfaceMesh_smooth = GM_smoothing(GMSurfaceMesh, smooth_value);
% Create mask of regions with no sensitivity
% And mask for plots (a little bit more strict)
[Jmask, Jmask2] = Jacobian_display_masking_speech(GMSurfaceMesh_smooth, J_GM_wav1);
% Mask nodes with no sensitivity
% use Jmask2 for plots improves boundary regions
Jpos = GMSurfaceMesh_smooth.node((Jmask==0),:);
Mpos = [Mpos; Jpos];
% Get index of channels in previous order and channels in hmr order

% Feed data into node positions
CorrV_HbO = zeros(length(Mpos),1);
CorrV_HbO(1:ch) = betas(1,:);
CorrV_HbR = zeros(length(Mpos),1);
CorrV_HbR(1:ch) = betas(2,:);
% Interpolation of values to the rest of the head
F_HbO = scatteredInterpolant(Mpos,CorrV_HbO,'linear', 'linear');
val_HbO = F_HbO(GMSurfaceMesh_smooth.node(:,1),GMSurfaceMesh_smooth.node(:,2),GMSurfaceMesh_smooth.node(:,3));
F_HbR = scatteredInterpolant(Mpos,CorrV_HbR,'linear', 'linear');
val_HbR = F_HbR(GMSurfaceMesh_smooth.node(:,1),GMSurfaceMesh_smooth.node(:,2),GMSurfaceMesh_smooth.node(:,3));
% Inverse distance weighted interpolation
% Looks good too, but requires tuning the parameters
%val_HbO = idw(Mpos, CorrV_HbO, GMSurfaceMesh.node, p, inf, L);
%val_HbR = idw(Mpos, CorrV_HbR, GMSurfaceMesh.node, p, inf, L);
% Define positions far from channels as zero, not plotting
% Improves visualization at the boundaries
d=20;
for i = 1:ch
    NaN_index(:,i) = Mpos(i,1)-d <= GMSurfaceMesh.node(:,1) &  GMSurfaceMesh.node(:,1) <= Mpos(i,1) + d &...
        Mpos(i,2)-d <= GMSurfaceMesh.node(:,2) &  GMSurfaceMesh.node(:,2) <= Mpos(i,2) + d &...
        Mpos(i,3)-d <= GMSurfaceMesh.node(:,3) &  GMSurfaceMesh.node(:,3) <= Mpos(i,3) + d ;
end
NaN_index = sum(NaN_index,2)>0;
val_HbO_final = val_HbO;
val_HbO_final(~NaN_index) = 0;
val_HbR_final = val_HbR;
val_HbR_final(~NaN_index) = 0;
% ADD THRESHOLDING plot 80% higher
% th_HbO = sort(abs(val_HbO_final), 'descend');
% th_HbO(th_HbO==0) = [];
% th_HbO = th_HbO(ceil(length(th_HbO)*0.5));
%
% th_HbR = sort(abs(val_HbR_final), 'descend');
% th_HbR(th_HbR==0) = [];
% th_HbR = th_HbR(ceil(length(th_HbR)*0.5));
%
% pos_th_HbO = (val_HbO_final> th_HbO).*val_HbO_final;
% neg_th_HbO = (val_HbO_final< -th_HbO).*val_HbO_final;
%
% pos_th_HbR = (val_HbR_final> th_HbR).*val_HbR_final;
% neg_th_HbR = (val_HbR_final< -th_HbR).*val_HbR_final;
%
% val_HbO_final = pos_th_HbO + neg_th_HbO;
% val_HbR_final = pos_th_HbR + neg_th_HbR;

% PLOT FIGURE
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
for i=1:2
    subplot(2,2,i);
    % Plot head
    hh = trisurf(GMSurfaceMesh.face,GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3),val_HbO_final.*Jmask2');hold on
    set(hh,'edgecolor','none','facecolor','interp','facelighting','phong')
    set(hh,'diffusestrength',0.9,'specularstrength',.1,'ambientstrength',.6);
    h = light;
    lighting gouraud
    %lightangle(h, 15, 100)       
    xlabel('x [mm]'),ylabel('y [mm]'),zlabel('z [mm]'); daspect([1 1 1])
    caxis([-0.15 0.15]);%colorbar
    % Plot channels
    plot3(Mpos(1:ch,1),Mpos(1:ch,2),Mpos(1:ch,3),'ko','markerfacecolor','k','markersize',4);
    plot3(Mpos(seed,1),Mpos(seed,2),Mpos(seed,3),'ro','markerfacecolor','r','markersize',6);
    % Display 3 different views
    if i == 1
        view([180 0 0]),axis tight,axis off    
        light('position', [1 0 0], 'style', 'infinite')
        title('HbO - LH ', 'fontsize', 20)
    else
        view([180*(-1) 0 0]),axis tight,axis off
        light('position', [-1 0 0], 'style', 'infinite')
        title('HbO - RH ', 'fontsize', 20)
    end
    colormap(greyJet)
end

hsup = suptitle(['PPI seed channel: ' num2str(seed)]);
set(hsup, 'fontsize',24)

for i=1:2
    subplot(2,2,i+2)
    % Plot head
    hh = trisurf(GMSurfaceMesh.face,GMSurfaceMesh.node(:,1),GMSurfaceMesh.node(:,2),GMSurfaceMesh.node(:,3),val_HbR_final.*Jmask2');hold on
    set(hh,'edgecolor','none','facecolor','interp','facelighting','phong')
    set(hh,'diffusestrength',0.9,'specularstrength',.1,'ambientstrength',.6);
    h = light;
    lighting gouraud
    %lightangle(h, 15, 100)    
    % Plot channels
    plot3(Mpos(1:ch,1),Mpos(1:ch,2),Mpos(1:ch,3),'ko','markerfacecolor','k','markersize',4);
    plot3(Mpos(seed,1),Mpos(seed,2),Mpos(seed,3),'ro','markerfacecolor','r','markersize',6);
    %Properties
    xlabel('x [mm]'),ylabel('y [mm]'),zlabel('z [mm]'),daspect([1 1 1])
    caxis([-0.15 0.15]);%colorbar
    % Display 3 different views
    if i == 1
        view([180 0 0]); axis tight; axis off
        light('position', [1 0 0], 'style', 'infinite')
        title('HbR - LH ', 'fontsize', 20)
        
    else
        view([180*(-1) 0 0]),axis tight,axis off
        light('position', [-1 0 0], 'style', 'infinite')
        colorbar
        set(gca, 'fontsize', 32)
        title('HbR - RH ', 'fontsize', 20)
        
    end
    colormap(greyJet)
end

figure
imagesc(betas, [-0.15 0.15]); colormap(greyJet)

% Save figure






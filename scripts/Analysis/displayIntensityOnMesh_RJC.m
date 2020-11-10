function [h l]=displayIntensityOnMesh_RJC(mesh, intensity)

% Usage:
%
% h=displayIntensityOnMesh(mesh,intensity,thresholds,axes_order)
%
% Written by Jay Dubb amended by RJC to maintain axes and orientation
%

% To display my results on this format I need
% - Age appropiate GM mesh.
% - The position of my channels on the GM mesh
% - group ICA coefficients or Z-scores (same size as nodes)
% the nodes associated with the channels will have a value or interpolated
% value, and the rest will be zero
nodes = mesh.node;
elem  = mesh.face;

% find nodes at distance less than 10mm from node_center
%dist = 10;
%node_center1 = nodes(27,:);
%roi1 = sqrt(sum((nodes - repmat(node_center1,length(nodes),1)).^2,2));
%nodesWithin1 = roi1<=dist;

% find homologous node based on coordinates (x opposite sign)
%node_center1_RH = node_center1; node_center1_RH(1) = -node_center1_RH(1);

% find closest node
%M = size(nodes,1);
%dists = sqrt(sum((nodes - repmat(node_center1_RH,M,1)).^2,2));
%[~ ,ind] = min(dists);
%node_center1_RH = nodes(ind,:);
%node_center2 = node_center1_RH;
%roi2 = sqrt(sum((nodes - repmat(node_center2,length(nodes),1)).^2,2));
%nodesWithin2 = roi2<=dist;

%intensity = zeros(1, length(nodes));
%intensity(nodesWithin1) = 1;
%intensity(nodesWithin2) = 1;


h=trisurf(elem(:,1:3), nodes(:,1), nodes(:,2), nodes(:,3),intensity,'facecolor','interp','edgealpha',0.1);
set(h,'diffusestrength',0.9,'specularstrength',.1,'ambientstrength',.35);
set(h,'Facelighting','phong');
light('Position',[0 0 1],'Style','infinite');
light('Position',[0 0 -1],'Style','infinite');
%light('Position',[0 1 1],'Style','infinite');
% light('Position',[1 1 0],'Style','infinite');
light('Position',[0 0 -1],'Style','infinite');
axis equal;axis off;

%Set balanced colorbar;
%colorbar; 
colormap jet;
climits = caxis;
limit = max(abs(climits));
caxis([-limit limit]);
%pause(0.5)



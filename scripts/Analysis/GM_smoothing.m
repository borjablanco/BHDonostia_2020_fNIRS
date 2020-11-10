function GMSurfaceMesh_smooth = GM_smoothing(GMSurfaceMesh, smooth_value)
GMSurfaceMesh_smooth = GMSurfaceMesh;

for i = 1:length(GMSurfaceMesh.node)
    
    idx = getNodesWithinDist_RJC(GMSurfaceMesh.node(i,:),GMSurfaceMesh.node, smooth_value);
    idx_nodes =   GMSurfaceMesh.node(idx,:);
    new_coord = mean(idx_nodes);
    
    GMSurfaceMesh_smooth.node(i,:)  = new_coord;
end

%figure
%subplot(121)
%displayIntensityOnMesh_RJC(GMSurfaceMesh_smooth, log10(JsumNorm));
%caxis([-3 0])
%subplot(122)
%displayIntensityOnMesh_RJC(GMSurfaceMesh, log10(JsumNorm));
%caxis([-3 0])
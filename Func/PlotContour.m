function PlotContour(node,elem,value, isMesh)
% % show the contour
% % ** code by P.M.H @bit.edu.cn (CN) **
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com
% %  ---------------------------------------
if nargin < 4
    isMesh = 0;
end
% figure
if isMesh
    patch('Faces',elem(:,1:4),'Vertices',node,'FaceVertexCData',value,'FaceColor','interp');
else
    patch('Faces',elem(:,1:4),'Vertices',node,'FaceVertexCData',value,'FaceColor','interp','EdgeColor','non');
end

colorbar
% caxis([0 1])

title({['max value: ',num2str(max(value))],['min value: ',num2str(min(value))]})

end
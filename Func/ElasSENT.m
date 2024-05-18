function BoundaryCondition = ElasSENT(fixNode, NDof, inc)
% SENT - Single Edge Notched Tension
% % ** code by P.M.Hu @bit.edu.cn (CN) **
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com

dsipBC = inc;
fixNode(fixNode(:,3)>0, 3) = dsipBC;
MoveNode = fixNode(fixNode(:,3)>0, 1:2);

% Neumann(Natural) Boundary Condition
F = zeros(NDof, 1);

% BCs
BoundaryCondition.DirchletDOF = fixNode(:,1)*2 + (fixNode(:,2)-2);
BoundaryCondition.Dirichlet   = fixNode(:,3);
BoundaryCondition.FreeDOF     = setdiff([1:NDof]', BoundaryCondition.DirchletDOF);
BoundaryCondition.BDforce     = MoveNode(:,1)*2 + (MoveNode(:,2)-2);
BoundaryCondition.RHS         = F;

end
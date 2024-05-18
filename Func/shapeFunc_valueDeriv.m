function [GaussInfo] = shapeFunc_valueDeriv(elem, node, Para)
% -------------------------------------------------------------------
% Precomputing the shape function's values and derivatives at Gauss Pts
% For: physical value at GPt (shape function's values * nodal values)
%      strain-disp matrix at GPt (shape function's derivatives)
% ---------------------------------------------------------------------
% %  Create by P.M.Hu @ bit.edu.cn
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com

GaussInfo = struct;
Dof = Para.ndim;

numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele

NGPs = [2, 2]; % full integration
[gp, wgt] = gauss_quadrature(NGPs(1), NGPs(2));

EleShapeDeriv = cell(numEle, 1);  % Ele: shape function's deriv at GPts
EleShapeVal = cell(numEle, 1);    % Ele: shape function's value at GPts
JW = cell(numEle, 1);             % Ele: det(Jacobian) * weight at GPts
Hisy = cell(numEle, 1);           % Ele: reference energy at GPts

for ei = 1 : numEle
    Ndcoord = node(elem(ei,:),:);
    
    dRdxGaussPt = zeros(Dof, numEleNd, size(gp,1)); 
    RGaussPt = zeros(size(gp,1), numEleNd);
    JacWeight = zeros(size(gp,1),1);
    
    HisyGaussPt = zeros(size(gp,1),1); % predefine crack with nodal value, HisyGaussPt=0
    
    for gpti = 1 : size(gp,1)
        GPtRef = gp(gpti,:);  % reference parametric coordinates for each integration point
        W = wgt(gpti);        % weigths for each integration point
        
        [R, R1] = ShapeDerivatives2D(GPtRef(1),GPtRef(2));
        dxdxi = R1 * Ndcoord; % Jacobi matrix
        J = det(dxdxi);  % jacobian
        % Compute derivatives of basis functions w.r.t physical coordinates
        dRdx = dxdxi^(-1) * R1;
        
        dRdxGaussPt( :, :, gpti) = dRdx; %
        RGaussPt(gpti, :) = R;           %
        JacWeight(gpti) = J * W;
         
    end
    
    %
    EleShapeDeriv{ei} = dRdxGaussPt;
    EleShapeVal{ei} = RGaussPt;
    JW{ei} = JacWeight;
    Hisy{ei} = HisyGaussPt; % all zeros
end

GaussInfo.SpDeriv = EleShapeDeriv;
GaussInfo.SpVal = EleShapeVal;
GaussInfo.JW = JW;
GaussInfo.Hisy = Hisy;

end

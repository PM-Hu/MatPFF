function [Disp] = assembleElasKK(GaussInfo, elem, Phi, Para, BC)
% -------------------------------------------------------------------
% % ** code by P.M.Hu @BIT (CN) **
% Calculate stiffness matrices
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com
% -------------------------------------------------------------------

E = Para.E; % elastic modulus
nu = Para.nu; % poisson's ratio
if Para.isStress == 1  % plane stress
    D = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
else  % plane strain
    D = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))];
end

numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
numEDofs = numEleNd * Para.ndim;

KVals = zeros(numEDofs^2, numEle); % store the stiff matrix
for ei = 1 : numEle
        elei = elem(ei,:);
        Ke = zeros(numEDofs); % element stiff-u
        
        % loading FEM information
        dRdxGaussPt = GaussInfo.SpDeriv{ei};
        RGaussPt = GaussInfo.SpVal{ei};  
        JW = GaussInfo.JW{ei};
        
        for gpti = 1 : size(dRdxGaussPt,3)
            
            %Compute derivatives of basis functions w.r.t physical coordinates
            dRdx = dRdxGaussPt( :, :, gpti);
            R = RGaussPt(gpti, :);
            
            % %      _                                             _
            % %     | N_{1, x} 0         ... N_{m, x} 0         ...|
            % %  B =| 0        N_{1, y}  ... 0        N_{m, y}  ...|
            % %     | N_{1, y} N_{1, x}  ... N_{m, y} N_{m, x}  ...|
            % %      -                                             -
            % % \sigma = [\sigma_{xx} \sigma_{yy} \sigma_{xy}]
            B = zeros(3, 2 * numEleNd);
            B(1, 1 : 2 : 2 * numEleNd) = dRdx(1, :);
            B(2, 2 : 2 : 2 * numEleNd) = dRdx(2, :);
            
            B(3, 1 : 2 : 2 * numEleNd) = dRdx(2, :);
            B(3, 2 : 2 : 2 * numEleNd) = dRdx(1, :);
            
            GPphi = R * Phi(elei);  % phi at GPt
            
            StiffDegrad = (1-GPphi)^2 + 1d-12;  % gb = (1-d)^2
                        
            % compute element stiffness at quadrature point
            Ke = Ke + B' * D * B * JW(gpti) * StiffDegrad;  % hybrid model 
            
        end
        KVals(:, ei) = Ke(:);
end

J = repmat(1:numEDofs, numEDofs, 1);
I = J';
El = elem';
Eldofs = [El(:)*2-1 El(:)*2]';
ElConn = reshape(Eldofs(:),numEleNd*2, numEle);
ElConn = ElConn';

ii = ElConn(:, I(:))';
jj = ElConn(:, J(:))';

K = sparse(ii(:), jj(:), KVals(:)); 
K = (K + K')/2;

% solve
Disp = zeros(Para.NNd*Para.ndim,1);
F = zeros(Para.NNd*Para.ndim,1);
F(BC.FreeDOF) = F(BC.FreeDOF) - K(BC.FreeDOF, BC.DirchletDOF) * BC.Dirichlet;
Disp(BC.FreeDOF) = K(BC.FreeDOF, BC.FreeDOF) \ F(BC.FreeDOF);
Disp(BC.DirchletDOF) = BC.Dirichlet; % enforce BC 
end
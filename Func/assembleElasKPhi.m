function [Phi] = assembleElasKPhi(GaussInfo, elem, Disp, Para)
% -------------------------------------------------------------------
% % ** code by P.M.Hu @bit.edu.cn **
% Calculate Phi stiffness matrices
% ------------------------------------------------------------------

Gc = Para.Gc;    % energy release rate
Len = Para.Len;  % length scale
lambda = Para.lambda;
mu = Para.mu; % Lame Constant
PFModel = Para.PFModel;

switch PFModel  % % =======AT1/2 model parameters========
    case 1 % AT2 --Bourdin 
        dalpha = @(x) 2*x;
        ddalpha = 2;
        Calpha = 2;
        a2 = 1;
    case 2 % AT1 -- Pham
        dalpha = @(x) 1;
        ddalpha = 0;
        Calpha = 8/3;
        a2 = 1;
end

numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
numEDofs = numEleNd * Para.ndim;

KphiVals = zeros((numEleNd)^2, numEle); % store the stiff matrix
FphiVals = zeros(numEleNd, numEle); % store the rhs vector

for ei = 1 : numEle
    elei = elem(ei,:);
    eleDOFs = reshape([2*elei-1; 2*elei], numEDofs,1);
    
    Ke = zeros(numEleNd); % element stiff-phi
    Fe = zeros(numEleNd,1); % element rhs-phi
    
    % loading FEM information
    dRdxGaussPt = GaussInfo.SpDeriv{ei};
    RGaussPt = GaussInfo.SpVal{ei};
    JW = GaussInfo.JW{ei};
    HisyGaussPt = GaussInfo.Hisy{ei};
    
    for gpti = 1 : size(dRdxGaussPt,3)
        
        %Compute derivatives of basis functions w.r.t physical coordinates
        dRdx = dRdxGaussPt( :, :, gpti);
        R = RGaussPt(gpti, :);
        
        B = zeros(3, 2 * numEleNd);
        B(1, 1 : 2 : 2 * numEleNd) = dRdx(1, :);
        B(2, 2 : 2 : 2 * numEleNd) = dRdx(2, :);
        
        B(3, 1 : 2 : 2 * numEleNd) = dRdx(2, :);
        B(3, 2 : 2 : 2 * numEleNd) = dRdx(1, :);
        
        Bd = zeros(2, numEleNd);
        Bd(1, :) = dRdx(1, :);
        Bd(2, :) = dRdx(2, :);  % hat{B}
        
        GPstrain = B * Disp(eleDOFs); % strain at GPt (tn)
        GPHisy = HistoryVariable(GPstrain, HisyGaussPt(gpti), lambda, mu); % max reference energy
        switch PFModel
            case 1 % AT2
                % derivatives of degradation function (g - degradation function)
                d2gdd2 = 2;         
                dgdd = -2;           
            case 2 % AT1
                % derivatives of degradation function (g - degradation function)
                d2gdd2 = 2; 
                dgdd = -2; 
                GPHisy = max(GPHisy, 3*Gc/(16*Len)); 
        end
        
        % compute element stiffness at quadrature point
        Ke = Ke + Gc*2*a2*Len/Calpha * (Bd' * Bd) * JW(gpti);  %
        Ke = Ke + (Gc*ddalpha/(Calpha*Len) + d2gdd2*GPHisy) * (R' * R) * JW(gpti); % 
        
        Fe = Fe + (dgdd*GPHisy) * R' * JW(gpti); % 
        Fe = Fe + Gc/(Calpha*Len) *( R'*dalpha(0) )* JW(gpti);
    end
    KphiVals(:, ei) = Ke(:);
    FphiVals(:, ei) = Fe(:);
end

% %
J = repmat(1 : numEleNd, numEleNd, 1);
I = J';
ElConn = elem;

ii = ElConn(:, I(:))';
jj = ElConn(:, J(:))';

Kphi = sparse(ii(:), jj(:), KphiVals(:)); % assemble Kphi
Kphi = (Kphi + Kphi')/2;

kk = reshape(ElConn',[],1);
Fphi = sparse(kk, 1, FphiVals(:));  % assemble rhs

% Solve the system by mldivide - direct method
Phi = Kphi \ -Fphi; %

end
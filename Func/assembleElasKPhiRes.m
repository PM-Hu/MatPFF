function [Kphi, Fphi] = assembleElasKPhiRes(GaussInfo, elem, Disp, Phi, Para)
% -------------------------------------------------------------------
% % ** code by P.M.H. @BIT (CN) **
% % Calculate Phi stiffness matrices
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com
% %  ---------------------------------------

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
    case 3 % PFCZM -- Wu - Linear soften law
%         dalpha = @(x) 2 - 2*x;
%         ddalpha = 2*(1-2);
%         Calpha = pi;
%         a2 = 1;
%         p = 2;
%         b1 = 4*Para.Lch/(Calpha*Len);
%         b2 = -0.5;
%         b3 = 0;
        % PFCZM -- Wu - Cornelissen soften law
        dalpha = @(x) 2 - 2*x;
        ddalpha = 2*(1-2);
        Calpha = pi;
        a2 = 1;
        p = 2;
        b1 = 4*Para.Lch/(Calpha*Len);
        b2 = 1.3868;
        b3 = 0.6567;
end

numEleNd  = size(elem, 2);  % 单元结点数
numEle = size(elem, 1); % 单元数
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
        
        GPphi = R * Phi(elei);       % phi at GPt
        dGPphi = dRdx * Phi(elei);   % phi-first-dirv at GPt        
        
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
        
        Bd = zeros(2, numEleNd);
        Bd(1, :) = dRdx(1, :);
        Bd(2, :) = dRdx(2, :);  % hat{B}
        
        GPstrain = B * Disp(eleDOFs);    % strain at GPt (tn)
        GPHisy = WuHistoryVariable(GPstrain, HisyGaussPt(gpti), Para); % max reference energy
        switch PFModel
            case 1 % AT2
                % derivatives of degradation function (g - degradation function)
                d2gdd2 = 2;          % gb = (1-d)^2
                dgdd = -2*(1-GPphi); % gb = (1-d)^2
            case 2 % AT1
                % derivatives of degradation function (g - degradation function)
                d2gdd2 = 2;          % gb = (1-d)^2
                dgdd = -2*(1-GPphi); % gb = (1-d)^2
                GPHisy = max(GPHisy, 3*Gc/(16*Len)); 
            case 3 % PFCZM
                % derivatives of degradation function
                fac1    =  (1-GPphi)^p; % PFCZM
                dfac1   = -p*(1 - GPphi)^(p - 1);
                ddfac1  =  p*(p - 1)*(1 - GPphi)^(p - 2);
                
                fac2    =  fac1   + b1*GPphi + b1*b2*GPphi^2 + b1*b2*b3*GPphi^3;
                dfac2   =  dfac1  + b1 + 2*b1*b2*GPphi + 3*b1*b2*b3*GPphi^2;
                ddfac2  =  ddfac1 + 2*b1*b2 + 6*b1*b2*b3*GPphi;
                
                dgdd  =  (dfac1*fac2  - fac1*dfac2)/(fac2^2);
                d2gdd2  =  ((ddfac1*fac2 - fac1*ddfac2)*fac2 - 2*(dfac1*fac2 - fac1*dfac2)*dfac2)/(fac2^3);
                GPHisy = max(GPHisy, Gc/(2*Para.Lch)); 
        end
        
        % compute element stiffness at quadrature point
        Ke = Ke + Gc*2*a2*Len/Calpha * (Bd' * Bd) * JW(gpti);  %
        Ke = Ke + (Gc*ddalpha/(Calpha*Len) + d2gdd2*GPHisy) * (R' * R) * JW(gpti); % %
        
        Fe = Fe + (dgdd*GPHisy) * R' * JW(gpti); % % 
        Fe = Fe + Gc/(Calpha*Len) *( R'*dalpha(GPphi) + 2*a2*Len^2*(Bd'*dGPphi) )* JW(gpti);
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

end

% % MIT License
% % 
% % Copyright (c) 2024 PM-Hu
% % 
% % Permission is hereby granted, free of charge, to any person obtaining a copy
% % of this software and associated documentation files (the "Software"), to deal
% % in the Software without restriction, including without limitation the rights
% % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % copies of the Software, and to permit persons to whom the Software is
% % furnished to do so, subject to the following conditions:
% % 
% % The above copyright notice and this permission notice shall be included in all
% % copies or substantial portions of the Software.
% % 
% % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% % SOFTWARE.
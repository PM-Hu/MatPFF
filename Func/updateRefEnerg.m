function [GaussInfo, bn] = updateRefEnerg(GaussInfo, elem, Disp, Para)
% -------------------------------------------------------------------
% % ** code by P.M.H @bit.edu.cn (CN) **
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com
% ------------------------------------------------------------------
E = Para.E; % elastic modulus
nu = Para.nu; % poisson's ratio
lambda = Para.lambda;
mu = Para.mu; % Lame Constant

if Para.isStress == 1  % plane stress
    D = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
else  % plane strain
    D = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))];
end

numEleNd  = size(elem, 2);  % 
numEle = size(elem, 1); % 
numEDofs = numEleNd * Para.ndim;

bnvals = zeros(numEDofs, numEle);
for ei = 1 : numEle
        elei = elem(ei,:);
        eleDOFs = reshape([2*elei-1; 2*elei], numEDofs,1);
        be = zeros(numEDofs, 1);

        % loading FEM information
        dRdxGaussPt = GaussInfo.SpDeriv{ei};
        HisyGaussPt = GaussInfo.Hisy{ei};
        JW = GaussInfo.JW{ei};
        
        for gpti = 1 : size(dRdxGaussPt,3)
            
            %Compute derivatives of basis functions w.r.t physical coordinates
            dRdx = dRdxGaussPt( :, :, gpti);
            
            B = zeros(3, 2 * numEleNd);
            B(1, 1 : 2 : 2 * numEleNd) = dRdx(1, :);
            B(2, 2 : 2 : 2 * numEleNd) = dRdx(2, :);
            
            B(3, 1 : 2 : 2 * numEleNd) = dRdx(2, :);
            B(3, 2 : 2 : 2 * numEleNd) = dRdx(1, :);
           
            GPstrain = B * Disp(eleDOFs);    % strain at GPt (tn)     
            GPHisy = HistoryVariable(GPstrain, HisyGaussPt(gpti), lambda, mu); % max reference energy
            HisyGaussPt(gpti) = GPHisy; % update reference energy in history
            
            be = be + B' * D * GPstrain * JW(gpti);
        end
        
        GaussInfo.Hisy{ei} = HisyGaussPt; % update
        bnvals(:, ei) = be(:);
end 

El = elem';
Eldofs = [El(:)*2-1 El(:)*2]';
ElConn = reshape(Eldofs(:),numEleNd*2, numEle);
ElConn = ElConn';

kk = reshape(ElConn',[],1);
bn = sparse(kk, 1, bnvals(:));  % assemble rhs

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

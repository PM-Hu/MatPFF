function  Phi = NewtonItPhaseField(GaussInfo, elem, Disp, Phi_1, Para)
% 
% *** Newton Raphson (NR) approximation for nonlinear problem ***
% ------------------ Now take a brief review --------------------------
% Given a differentiable function F(x), the root x* (that is F(x*) = 0)
% is approximated by the incremental {\Delta x} and initial guess {x0}: 
% x* = x0 + \Delta x,  where \Delta x = - F(x0) / F'(x0).
% ---------------------------------------------------------------------
% Create by P.M.Hu @ BIT(CN) 2022-11-20
% Please feel free to contact us with any questions! 
% - Email: pm_hu@outlook.com
% %  ---------------------------------------

NRresidual = 1; 
tol = 1d-4;
Phi = Phi_1;
while NRresidual > tol
    
    [Kphi, Rphi] = assembleElasKPhiRes(GaussInfo, elem, Disp, Phi, Para);
    
    f =  -Rphi;
    % Solve the system by mldivide - direct method
    dphi = Kphi \ f; % 
 
    Phi = Phi + dphi;
    NRresidual = norm(dphi);
end

Phi(Phi >= 1) = 1; % enforce boundary constraint [0 1]
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
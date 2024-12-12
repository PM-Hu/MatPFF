function Hisy = WuHistoryVariable(strain, Hisold, Para)
% % ** code by P.M.H. @BIT (CN) **
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com
% %  ---------------------------------------
% % Last update: 2024-12-12

E = Para.E; % elastic modulus
nu = Para.nu; % poisson's ratio

if Para.isStress == 1  % plane stress
    D = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
else  % plane strain
    D = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))];
end

Stress = D * strain;

E1 = Stress(1);  % components
E2 = Stress(2);  % 
E12 = Stress(3); % 

M=sqrt( (E1-E2)^2 + 4 *E12^2 );
epsilon1 = 0.5*(E1+E2) + 0.5*M;  % eigen value
epsilon2 = 0.5*(E1+E2) - 0.5*M;  % 

MajorStress = max(epsilon1, epsilon2);
MajorStress = max(MajorStress, 0);

referEnerg = (MajorStress^2) / (2*Para.E);

Hisy = max(referEnerg, Hisold);
end

function val = H(x)

if x >= 0
    val = 1;
else
    val = 0;
end

end

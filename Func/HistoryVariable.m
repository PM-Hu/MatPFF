function Hisy = HistoryVariable(strain, Hisold, lambda, mu)
% %  ** code by P.M.Hu @BIT (CN) **
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com

E1 = strain(1);  % components
E2 = strain(2);  % 
E12 = strain(3)/2; % 

M=sqrt( (E1-E2)^2 + 4 *E12^2 );
epsilon1 = 0.5*(E1+E2) + 0.5*M;  % eigen value
epsilon2 = 0.5*(E1+E2) - 0.5*M;  % 

trE = E1 + E2;  % volume strain

referEnerg = 0.5*lambda*( H(trE) * trE )^2 + mu*( H(epsilon1)*epsilon1^2+ H(epsilon2)* epsilon2^2);

Hisy = max(referEnerg, Hisold);
end

function val = H(x)

if x >= 0
    val = 1;
else
    val = 0;
end

end
function [R, R1] = ShapeDerivatives2D(xi,eta)
% R - 形函数在高斯点处的值
% R1 - 形函数在高斯点处的导数

% N1 = (1-xi)(1-eta)/4
% N2 = (1+xi)(1-eta)/4
% N3 = (1+xi)(1+eta)/4
% N4 = (1-xi)(1+eta)/4

R(1) = (1-xi)*(1-eta)/4;
R(2) = (1+xi)*(1-eta)/4;
R(3) = (1+xi)*(1+eta)/4;
R(4) = (1-xi)*(1+eta)/4;

dNdxi = zeros(4,1);
dNdeta = zeros(4,1);

dNdxi(1) = (-1+eta)/4;
dNdxi(2) = (1-eta)/4;
dNdxi(3) = (1+eta)/4;
dNdxi(4) = (-1-eta)/4;

dNdeta(1) = (-1+xi)/4;
dNdeta(2) = (-1-xi)/4;
dNdeta(3) = (1+xi)/4;
dNdeta(4) = (1-xi)/4;

R1 = [dNdxi,dNdeta]';
end
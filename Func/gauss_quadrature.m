function [Q, W] = gauss_quadrature(noGpsX, noGpsY)
%  Calculate gauss integration points and its weights

a = -1;    % lower limit of integral interval
b = 1;     % upper limit of integral interval 

Q = zeros(noGpsX*noGpsY,2);
W = zeros(noGpsX*noGpsY,1);

[Q1,W1] = gauss_legendre(a,b,noGpsX);
[Q2,W2] = gauss_legendre(a,b,noGpsY);
n = 0;
for j = 1:noGpsY
    for i = 1:noGpsX
        n = n+1;
        Q(n,:) = [Q1(i),Q2(j)];
        W(n) = W1(i)*W2(j);
    end
end

end


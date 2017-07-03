function [ R ] = K2R( K,D,N,c,f_number,k,b )
%Radius R for an eye of certain volume  
%
% K volume divided by 4 pi /3
% D facet diameter
% N number of transduction units per ph
% c factor
% b multiplier for tr units in ph. 3 for Comp, 6 for CompN

LoadCalliphoraV0;
n3=1.35;

K1 = K;
L=N*k+D*f_number*n3;

%R=(3*L*L + sqrt(9*L*L*L*L - 12*(L+3*(4*pi/3)*(K0)*c*N/D/D/(2*sqrt(3)/3))*(L*L*L - K1)))/(6*(L+3*(4*pi/3)*(K0)*c*N/D/D/(2*sqrt(3)/3)));

%K = 3*L.*R.*R - 3*L.*L.*R + L.*L.*L  + K0*b*c*4*pi*N./D./D./(2*sqrt(3)/3).*R.*R;

A=3*L + K0*b*c*4*pi*N./D./D./(2*sqrt(3)/3);
B=-3*L.*L;
C=L.*L.*L - K1;

R=(-B+sqrt(B*B-4*A*C))/2/A;

end


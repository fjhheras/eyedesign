function [ K ] = R2K( R,D,N,c,f_number,k,b)
%Volume divided by pi3/4 of an eye with radious R
%
% 
% D facet diameter
% N number of transduction units per ph
% c factor
% b multiplier for tr units in ph. 3 for Comp, 6 for CompN

n3=1.35;

L=N*k+D*f_number*n3;

K = 3*L.*R.*R - 3*L.*L.*R + L.*L.*L + V_equiv_energy_cost(c,b,D,N,R); %K0*b*c*4*pi*N./D./D./(2*sqrt(3)/3).*R.*R;
    
  

end


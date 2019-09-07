function [ K ] = R2K( R,D,N,c,f_number,k,b,n)
%Volume divided by pi3/4 of an eye with radious R
%
% 
% D facet diameter
% N number of transduction units per ph
% c factor
% b multiplier for tr units in ph. 3 for Comp, 6 for CompN

L=N*k+cone_length(D,f_number,n);
K = 3*L.*R.*R - 3*L.*L.*R + L.*L.*L + V_equiv_energy_cost(c,b,D,N,R)./(4*pi/3) ; %K0*b*c*4*pi*N./D./D./(2*sqrt(3)/3).*R.*R;
    
  

end


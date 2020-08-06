function [ R ] = R_in_simple( V,L,d,c,N, b)
%R_IN_SIMPLE Calculate R in a simple eye, given
%   V - Volume equivalent of investment
%   L - Ph length
%   d - Distance between ph centers (normally Ph width)
%   c - cost
%   N - Number tr units per ph
%   b - factor 

  K1 = V / (4*pi/3);

  syms r; 
  
  %N_ph = b*4*pi*r*r/d/d/(2*sqrt(3)/3);
  %V_extra = c * N * N_ph; % Volume equivalent of transduction cost
  
  % V_equiv_energy_cost is proportional to R^2, so we can evaluate
  % in 1 and then multiply by r
  V_extra = V_equiv_energy_cost(c,b,d,N,r) %* r * r
  
  R = vpasolve( K1 - (r+L)^3 - V_extra / (4*pi/3) );
  
  
  %And we keep the only positive real root
  R = eval(R);
  R = R(imag(R)==0);
  R = R(R>0);
  
  if isempty(R)
      % If none, it means it is impossible to solve
      R = -1 
  end

end


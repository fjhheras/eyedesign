function [ R ] = R_in_simple( V,L,d,c,N,V0 )
%R_IN_SIMPLE Calculate R in a simple eye, given
%   V - Volume equivalent of investment
%   L - Ph length
%   d - Distance between ph centers (normally Ph width)
%   c - cost
%   N - Number tr units per ph

  K1 = V / (4*pi/3);
  %old_inv = V / V0;
  %n3=1.35;
  syms r;
  % Num photoreceptors: 4*pi*r*r/d/d/(2*sqrt(3)/3)
  % Investment (K1) must be equal to 
  R = vpasolve( K1 - (r+L)^3 - c*N*(4*pi*r*r/d/d/(2*sqrt(3)/3)) / (4*pi/3) );
  %And we keep the only positive real root
  R = eval(R);
  R = R(imag(R)==0);
  R = R(R>0);


end


function [ R ] = R_in_simple( V,L,d,c,N,V0 )
%R_IN_SIMPLE Summary of this function goes here
%   Detailed explanation goes here

  K1 = V / (4*pi/3);
  old_inv = V / V0;
  %n3=1.35;
  %And we keep the only positive real root
  syms r;
  R=vpasolve(K1-(r+L)^3-K1/old_inv*c*4*pi/d/d/(2*sqrt(3)/3)*r*r*N);
  R = eval(R);
  R = R(imag(R)==0);
  R = R(R>0);


end


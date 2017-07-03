function [ V_eq ] = V_equiv_energy_cost( c,b,D,N,R )
%V_EQUIV_ENERGY_COST Eq volume of the energy cost
%   D is lens diameter in compound eye, but rhabdom diameter in simple eye
LoadCalliphoraV0;
V_eq = V0*b*c*4*pi*N./D./D./(2*sqrt(3)/3).*R.*R;
          %CompN: 6*c*4*pi*y(2)/y(1)/y(1)/(2*sqrt(3)/3)*R*R
          %Simple:  c*4*pi*y(1)/d/d/(2*sqrt(3)/3)*R*R
end


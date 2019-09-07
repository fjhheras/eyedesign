function [ V_eq ] = V_equiv_energy_cost( c,b,D,N,R )
%V_EQUIV_ENERGY_COST Eq volume of the energy cost
%   D is lens diameter in compound eye, but rhabdom diameter in simple eye
% LoadCalliphoraV0;
%%Number of ommatidia is area of shpere divided by area of hex lenses
Nomm = 4*pi*R.*R./D./D./(2*sqrt(3)/3);
%%Number of transduction units is number of ommatidia times number of
%rabdomeres per ommatidia times number of t.u. per rhamdomere
Ntu = b*Nomm*N;
V_eq = c*Ntu;
          %CompN: 6*c*4*pi*y(2)/y(1)/y(1)/(2*sqrt(3)/3)*R*R
          %Simple:  c*4*pi*y(1)/d/d/(2*sqrt(3)/3)*R*R
end


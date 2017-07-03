S = NaturalImages(0.4,1)
LoadConstants
LoadCalliphoraV0

Da=[];
Na=[];
Ra=[];
Ha=[];
Hmaxa = [];
inv0a=[];
inv2as=[];
inva=logspace(-1,1,20)*V0;

%Nn=1e-4;
%c=power(10,-9.7);
c=0;

for inv =  inva
  K1 = inv/(4*pi/3)
  options = optimset('TolX',.0001);
  fH = @(x) -InfSimp( x, S, inv, c);  
  [y,H] = fminbnd(fH,0,nthroot(K1,3)/k,options);
  Na=[Na,y];
    
  L=y(1)*k;
  R = R_in_simple(inv,L,d,c,y(1),V0)
  old_inv = inv/V0
  inv0a = [inv0a,old_inv - V_equiv_energy_cost(c,1,d,y(1),R) %c*4*pi/d/d/(2*sqrt(3)/3)*R*R*y(1)];
  
  f = R/n3;
  D = f/f_number;
  Deltaphi = d/f; %rad
  
  Ra=[Ra, R];
  Da=[Da, D];
  Ha=[Ha, -H];
  inv2as=[inv2as (R+L)*(R+L)*(R+L)/K0]
  
end

La = Na*k

%N = 3000;
%Nt = 250;
%[S, fx, ft, Deltax, Deltat] = NatImagbigSH(N,Nt,0.4,1); %give me natural images stats with contrast 0.4 and speed 3 rad/s
S = NaturalImages(0.4,1)
%Narray = logspace(1,5,20);

Da=[];
Na=[];
Ra=[];
Ha=[];
Hmaxa = [];
inv0a=[];
inv2as=[];

%Nn=1e-4;
%c=power(10,-9.7);
%c=0;
%inv = 1;
%for c=logspace(-10.5,-8.5,20)
for inv = logspace(-1,1,20)  
   R0 = 31/(2*pi/180);
   L0 = 250 + 60; %um
   K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0; %it's a way of calculating the same material than in a compound
   K0 = K0*inv/2;
   k=4.3e-3;

  options = optimset('TolX',.0001);
  fH = @(x) -InfSimp( x, S, inv, c);  
  [y,H] = fminbnd(fH,0,nthroot(2*K0,3)/k,options);
  Na=[Na,y];
    
  f_number = 60/31;
  d=1.9;
  R0 = 31/(2*pi/180);
  L0 = 250 + 60; %um
  K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0; %it's a way of calculating the same material than in a compound eye
  K1 = K0*inv
  k = 4.3e-3
  L=y(1)*k;
  n3=1.35;
    
  %R = nthroot(K0,3)-L;
  
  syms r;
  R=vpasolve(K1-(r+L)^3-K1/inv*c*4*pi/d/d/(2*sqrt(3)/3)*r*r*y(1))
  R = eval(R)
  R = R(imag(R)==0);
  R = R(R>0)
  
  inv0a = [inv0a,inv - c*4*pi/d/d/(2*sqrt(3)/3)*R*R*y(1)];
  
  f = R/n3;
  D = f/f_number;
  Deltaphi = d/f; %rad
  
  Ra=[Ra, R];
  Da=[Da, D];
  Ha=[Ha, -H];
  inv2as=[inv2as (R+L)*(R+L)*(R+L)/K0]
  %Hmaxa = [Hmaxa, Hmax]
  
end

La = Na*k

%pause
% options = optimset('TolX',.00001);
% 
% fH = @(x) -InfNTu31( x(1), x(2), N, Deltax, Deltay, Deltat, fx,fy,ft, S, inv );
% y = fminsearch(fH,[Dmax, Nmax], options);
% Da=[Da,y(1)];
% Na=[Na,y(2)];
% 
% f_number = 31/60;
% R0 = 31/(3*pi/180);
% L0 = 285 + 4; %um
% K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0;
% k = 4e-3
% K0 = K0*inv
% L=y(2)*k+y(1)/f_number;
% R=(3*L*L + sqrt(9*L*L*L*L - 12*L*(L*L*L - K0)))/(6*L)
% 
% Ra=[Ra,R];
% 
% end
% 
% pa=Da.*Da./Ra
% inv=0.5:0.05:1.4

%surf(Narray, Darray,H./cost.*Deltaphi.*Deltaphi)
%surf(Narray, Darray,H) %if we fix cost by fixing volume, dividing by cost again is not necessary (and it's wrong cos R is no longer ctant)
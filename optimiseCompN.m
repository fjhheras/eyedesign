%clear all;
Sf = NaturalImages(0.4,1)
LoadConstants
LoadCalliphoraV0

DPoints=10;
NPoints=40;

Da=[];
Na=[];
Ra=[];
Ha=[];
Hmaxa = [];
inv0a = [];
inv2an=[];

Darray = linspace(15,40,DPoints);
Narray = logspace(4,5.5,NPoints);
cagada=0;
%Narray2 = linspace(3,6,50)
Harray2 = zeros(20,50)
inva = logspace(-1,1,20) * V0

invi=0;

%c=0.7e-10
%c=1e-10
%c=10e-10
c=0.0

for inv = inva
invi=invi+1;
Hmax=0;
inv
Dmax =-1;
Nmax =-1;
for iD=1:DPoints
    for iN = 1:NPoints
        H(iD,iN) = InfCompN(Darray(iD),Narray(iN), Sf, inv, c);
        if H(iD,iN)>Hmax
            Hmax = H(iD,iN);
            Dmax=Darray(iD);
            Nmax=Narray(iN);
        end
    end
end

options = optimset('TolX',.0001);
fH = @(x) -InfCompN( x(1), x(2),Sf, inv,c );
[y,H] = fminsearch(fH,[Dmax, Nmax], options);


if (InfCompN( Dmax, Nmax, Sf, inv,c ) > InfCompN( y(1), y(2),Sf, inv, c))
    cagada = cagada+1
    y(1)=Dmax;
    y(2)=Nmax;
end

Da=[Da,y(1)];
Na=[Na,y(2)];

%R0 = 31/(2*pi/180);
%L0 = 250 + 60; %um
%K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0;
%K1 = K0*inv

R = K2R(inv/(4*pi/3),y(1),y(2),c,f_number,k,6);
inv0= inv - V_equiv_energy_cost(c,6,y(1),y(2),R ); %6*c*4*pi*y(2)/y(1)/y(1)/(2*sqrt(3)/3)*R*R
inv2 = R2K(R,y(1),y(2),0,f_number,k,6)/K0; %(R*R*R-(R-L)*(R-L)*(R-L))/K0


Ra=[Ra,R];
Ha=[Ha,-H];
Hmaxa = [Hmaxa, Hmax];
inv0a=[inv0a,inv0];
inv2an=[inv2an,inv2];
end


%c=logspace(-7,-5,40)
pa=Da.*Da./Ra;
%inv=0.1:0.2:3
n3=1.35;
La = Na*k + Da*f_number*n3;

%surf(Narray, Darray,H./cost.*Deltaphi.*Deltaphi)
%surf(Narray, Darray,H) %if we fix cost by fixing volume, dividing by cost again is not necessary (and it's wrong cos R is no longer ctant)

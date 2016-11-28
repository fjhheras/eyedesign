%clear all;

DPoints=200;
NPoints=300;

Da=[];
Na=[];
Ra=[];
Ha=[];
Hmaxa = [];
inv0a=[];
inv2ac=[];

%N = 3000;
%Nt = 250;
%[S, fx, ft, Deltax, Deltat] = NatImagbigSH(N,Nt,0.4,3); %give me natural images stats with contrast 0.4 and speed 3 rad/s

Darray = linspace(0.1,20,DPoints);
Narray = logspace(1,5,NPoints);

Narray2 = linspace(3,6,50)
Harray2 = zeros(20,50)
%Nn=1e-4;
%inv = 1;
%c=power(10,-9.7);
%c=0;
invi=0;
%for c=logspace(-10.5,-8.5,20)
for inv=logspace(-4.6,0,20)
invi=invi+1;
Hmax=0
Dmax =-1;
Nmax =-1;
for iD=1:DPoints
    iD
    for iN = 1:NPoints
        H(iD,iN) = InfComp(Darray(iD),Narray(iN), N, Nt, Deltax, Deltat, fx, ft, S, inv, k,f_number,c, Nn);
        if H(iD,iN)>Hmax
            Hmax = H(iD,iN);
            Dmax=Darray(iD);
            Nmax=Narray(iN);
        end
    end
end


options = optimset('TolX',.0001);
fH = @(x) -InfComp( x(1), x(2), N, Nt, Deltax, Deltat, fx,ft,S, inv,k,f_number,c, Nn );
[y,H] = fminsearch(fH,[Dmax, Nmax], options);


if (InfComp( Dmax, Nmax, N, Nt, Deltax, Deltat, fx,ft,S, inv, k,f_number,c, Nn ) > InfComp( y(1), y(2), N, Nt, Deltax, Deltat, fx,ft,S, inv, k,f_number,c, Nn ))
    cagada = 1
    pause
end

Da=[Da,y(1)];
Na=[Na,y(2)];

%f_number = 60/31;
R0 = 31/(2*pi/180);
L0 = 250 + 60; %um
K0 = 3*L0*R0*R0 - 3*L0*L0*R0 + L0*L0*L0;
%k = 4.3e-3
K1 = K0*inv
%L=y(2)*k+y(1)*f_number;
%R=(3*L*L + sqrt(9*L*L*L*L - 12*(L+3*(4*pi/3)*(K0)*c*y(2)/y(1)/y(1)/(2*sqrt(3)/3))*(L*L*L - K1)))/(6*(L+3*(4*pi/3)*(K0)*c*y(2)/y(1)/y(1)/(2*sqrt(3)/3)));

R=K2R(K1,y(1),y(2),c,f_number,k,3);

inv0=inv-3*c*4*pi*y(2)/y(1)/y(1)/(2*sqrt(3)/3)*R*R;
inv2=R2K(R,y(1),y(2),0,f_number,k,3)/K0; %(R*R*R-(R-L)*(R-L)*(R-L))/K0

% for Ni = 1:50
%     
%     Harray2(invi,Ni) = InfCostHexBigSHlin( y(1), Narray2(Ni), N, Nt, Deltax, Deltat, fx,ft,S, inv, c, Nn ) / (-H);
%     
% end

Ra=[Ra,R];
Ha=[Ha,-H];
Hmaxa = [Hmaxa, Hmax];
inv0a=[inv0a inv0];
inv2ac=[inv2ac inv2];
end


%c=logspace(-7,-5,40)
pa=Da.*Da./Ra;
%inv=0.1:0.2:3
n3=1.35;
La = Na*k + Da*f_number*n3;

%surf(Narray, Darray,H./cost.*Deltaphi.*Deltaphi)
%surf(Narray, Darray,H) %if we fix cost by fixing volume, dividing by cost again is not necessary (and it's wrong cos R is no longer ctant)
function H = InfSimp( Ntu, N, Nt, Deltax, Deltat,fx,ft,S,inv, c, Nn)
%Information in a simple eye optimised for biglens
%Temporal resolution of photoreceptor not included!!!
%D = 31 %um
f_number = 60/31;
d = 1.9; %um
lambda = 0.5; %wave-lenght of yellow light in um

R0 = 31/(2*pi/180);

L0 = 250 + 60; %um
K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0; %it's a way of calculating the same material than in a compound
K0 = K0*inv; %half eye

Deltaphi0=31/R0;

n3=1.35;

k = 4.3e-3;
L= Ntu*k;

if ((L*L*L > K0)|(L<0))
    %L
    %K0
    H = 0;
    Deltaphi = 1;
    %pause
else
    
%R = nthroot(K0,3)-L;

syms R;
R=solve(K0-(R+L)^3-K0/inv*c*4*pi/d/d/(2*sqrt(3)/3)*R*R*Ntu)
eval(R)
%
R=ans(1);
f = R/n3;
D = f/f_number;
Deltaphi = d/f %rad 

H=InfG2( D, Ntu, N, Nt, Deltax, Deltat,fx,ft,S,R,f,Deltaphi, c, Nn, 0,1)

%H = H/(inv+c*4*pi*Ntu/Deltaphi/Deltaphi/(2*sqrt(3)/3));
H=H*2*sqrt(3)/pi/pi;
%StoN=sqrt(Tsignal/Tnoise)
end
end

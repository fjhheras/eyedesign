function H = InfComp (D, Ntu, N, Nt, Deltax, Deltat,fx,ft, S, inv,k,f_number, c, Nn)
%Information per cost considering frequencies in first Brioullin zone and keeping a
%constant volume modified by inv. Cost are modified by a parameter ,
%giving some cualitative factor for the energy costs of the photoreceptors
%D = 31 %um
%f_number = 60/31;
d = 1.9; %um
lambda = 0.5; %wave-lenght of yellow light in um
Deltaphi0 = 2*pi/180;
R0 = 31/(Deltaphi0);
L0 = 250 + 60; %um
K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0;
K1 = K0*inv
n3=1.35;


%k = 4.3e-3

L=Ntu*k+D*f_number*n3;


if ((K2Rdelta(K1,D,Ntu,c,f_number,k,3) < 0)|(L<0))
    H = 1;
    Deltaphi = 1;
else
    
R=K2R(K1,D,Ntu,c,f_number,k,3); %R=(3*L*L + sqrt(9*L*L*L*L - 12*(L+3*(4*pi/3)*(K1)*c*Ntu/D/D/(2*sqrt(3)/3))*(L*L*L - K0)))/(6*(L+3*(4*pi/3)*(K1)*c*Ntu/D/D/(2*sqrt(3)/3)))

if (R < L+R*2*d/D) %contour condition
    H = 0;
    Deltaphi = 1;
else

Deltaphi = D/R; %rad
f=D*f_number;

H=InfG2( D, Ntu, N, Nt, Deltax, Deltat,fx,ft,S,R,f,Deltaphi, c, Nn, 1,3)



%%H/(4*pi/Deltaphi/Deltaphi/(3*sqrt(3)/8))
%H=H/(inv+3*c*4*pi*Ntu/Deltaphi/Deltaphi/(2*sqrt(3)/3));
%%H = trapz(Ha)/(pi*K0+c*6*4*pi/Deltaphi/Deltaphi*k*Ntu*pi*d*d/4)
H=H*2*sqrt(3)/pi/pi;
end
end
end


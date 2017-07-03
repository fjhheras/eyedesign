function H = InfSimp( Ntu, S_f, inv, c)
%Information in a simple eye optimised for biglens
%Temporal resolution of photoreceptor not included!!!

LoadConstants
LoadCalliphoraV0

%f_number = 60/31;
%d = 1.9; %um
%lambda = 0.5; %wave-lenght of yellow light in um
%R0 = 31/(2*pi/180);
%L0 = 250 + 60; %um
%K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0; %it's a way of calculating the same material than in a compound
%K1 = K0*inv; %half eye
%Deltaphi0=31/R0;

K1 = inv/(4*pi/3)
old_inv = inv/V0
L= Ntu*k;

if ((L*L*L > K1)|(L<0))
    %L
    %K0
    H = 0;
    Deltaphi = 1;
    %pause
else
    
%R = nthroot(K0,3)-L;
R = R_in_simple( inv,L,d,c,Ntu,V0 )

%syms r;
%R=vpasolve(K1 -(r+L)^3 - K1/old_inv* c *4*pi/d/d/(2*sqrt(3)/3)*r*r*Ntu);
%And we keep the only positive real root
%R = eval(R);
%R = R(imag(R)==0);
%R = R(R>0);
%
f = R/n3
D = f/f_number
Deltaphi = d/f %rad 

H=InformationRate( D, Ntu, S_f,f,Deltaphi, 1)
%H=H*2*sqrt(3)/pi/pi;
end
end

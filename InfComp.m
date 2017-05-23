function H = InfComp (D, Ntu, S, inv, c)
%Information per cost considering frequencies in first Brioullin zone and keeping a
%constant volume modified by inv. Cost are modified by a parameter ,
%giving some cualitative factor for the energy costs of the photoreceptors

LoadConstants

Deltaphi0 = 2*pi/180;
R0 = 31/(Deltaphi0);
L0 = 250 + 60; %um
K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0;
K1 = K0*inv;

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

        H = InformationRate(D, 3*Ntu, S,f,Deltaphi, 1);
        %H=H*2*sqrt(3)/pi/pi;
    end
end
end


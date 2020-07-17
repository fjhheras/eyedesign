function H = InfComp (D, Ntu, S, inv, c)
% Output: bits per Hz and per sr

LoadConstants


K1 = inv/(4*pi/3) % It is easier to work with this number, rather than with volume
%L=Ntu*k+D*f_number*n3;
L=Ntu*k+cone_length(D,f_number,n3);


if ((K2Rdelta(K1, D, Ntu, c, f_number, k, 3) < 0)|(L<0))
    H = -1;
else
    
    R=K2R(K1,D,Ntu,c,f_number,k,3);
    
    if (R < L+R*2*d/D) %contour condition
        H = -1;
    else
        Deltaphi = D/R; %rad
        f=D*f_number;
        H = InformationRate(D, 3*Ntu, S,f,Deltaphi, 1);
    end
end
end


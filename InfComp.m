function H = InfComp (D, Ntu, S, inv, c, eye_type)
% Output: bits per Hz and per sr

LoadConstants

if eye_type == 'fly'
    b = 6 % Num rhabdomere in rhabdom
    contour_mult = 3
elseif eye_type == 'app'
    b = 3
    contour_mult = 1
else
    throw(MException('Wrong eye_type %s', eye_type))
end
    
K1 = inv/(4*pi/3) % It is easier to work with this number, rather than with volume
L=Ntu*k+cone_length(D,f_number,n3);


if ((K2Rdelta(K1, D, Ntu, c, f_number, k, b, n3) < 0)|(L<0))
    H = -1;
else
    
    R=K2R(K1, D, Ntu, c, f_number, k, b, n3);
    
    if (R < L+R*contour_mult*d/D) %contour condition
        H = -1;
    else
        Deltaphi = D/R; %rad
        f=D*f_number;
        H = InformationRate(D, b*Ntu, S,f,Deltaphi, 1);
    end
end
end


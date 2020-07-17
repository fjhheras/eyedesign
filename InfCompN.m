function H = InfCompN (D, Ntu, S, inv, c)
% Output: bits per Hz and per sr


LoadConstants

K1 = inv/(4*pi/3); % It is easier to work with this number, rather than with volume
L=Ntu*k+cone_length(D,f_number,n3);


if ((K2Rdelta(K1, D, Ntu, c, f_number, k, 6, n3) < 0)|(L < 0)) %3,6
    % Invalid geometrical configuration
    H = -1;
else
    
    R=K2R(K1,D,Ntu,c,f_number,k,6,n3); 
    
    if (R < L+R*3*d/D)
        H = -1;
        % Invalid geometrical configuration
    else
        Deltaphi = D/R; %rad
        f=D*f_number;
        H = InformationRate(D, 6*Ntu, S,f,Deltaphi, 1);
    end
end
end


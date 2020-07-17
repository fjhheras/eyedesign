function H = InfSimp( Ntu, S_f, inv, c)
%Information in a simple eye
% Output: bits per Hz and per sr

LoadConstants

K1 = inv/(4*pi/3)
%There is no cone to include in length for this case
L= Ntu*k;

if ((L*L*L > K1)|(L<0))
    H = -1;
    Deltaphi = 1;
    %pause
else

    R = R_in_simple(inv, L, d, c, Ntu, V0)
    f = R/n3 % focal length
    D = f/f_number % Lens diameter from f-number and focal length
    Deltaphi = d/f %rad 

    H=InformationRate( D, Ntu, S_f,f,Deltaphi, 1)

end
end

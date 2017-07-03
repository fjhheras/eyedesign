function H = snyderH( R, D, F_number, C , N )
%SNYDERH J. Howard and A. W. Snyder "Transduction as a limitation on compound eye
% function and design" (1983), eq 24, in the limit of large Smax and with the 
% substitution N = I*D^2*delta_phi^2 
%   R - radius
%   D - lens diameter
%   C - contrast of natural images
%   N - light intensity in ph/s

lambda = 0.5 %wave-length of yellow light in um
q = 0.5 %quantum capture efficiency, not in H and S, but reasonable
d=1.9
t_shutter = 1/71/2; %71Hz 
constant_info_array = 2*sqrt(3) %In HS (1983) they use 4 because they use sq array

f=F_number*D;
delta_phi=D/R;
ang_rab = d/f;
D=delta_phi*R;
ang_acc = sqrt(ang_rab.^2 + (lambda/D).^2);

M = @(mu) exp(-3.56*mu.^2*ang_acc.^2);
S2 = @(mu) (q*N*t_shutter*C.^2*(M(mu)).^2)/2;
F = @(mu) mu.*log2(1+S2(mu));

H = constant_info_array/t_shutter * quad(F,0,1/delta_phi/sqrt(3));
end


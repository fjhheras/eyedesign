function H = snyderH( R, D, C , I )
%SNYDERH J. Howard and A. W. Snyder "Transduction as a limitation on compound eye
% function and design" (1983) 
%   R - radius
%   D - lens diameter
%   C - contrast of natural images
%   I - light intensity

lambda = 0.5 %wave-length of yellow light in um
q = 0.5 %quantum capture efficiency
Smax = 1e5 %100 maximum signal-to-noise ratio
d=1.9
t_shutter = 1/71/2; 
f=60*D/31;
delta_phi=D/R;
ang_rab = d/f;
D=delta_phi*R;
ang_acc = sqrt(ang_rab.^2 + (lambda/D).^2);

M = @(mu) exp(-3.56*mu.^2*ang_acc.^2);
S2 = @(mu) (I*D^2*ang_rab.^2*t_shutter*Smax.^2*C.^2*(M(mu)).^2)/2/(I*D.^2*ang_rab.^2*t_shutter + Smax.^2);
F = @(mu) mu.*log(1+S2(mu));

H = 4/t_shutter*quad(F,0,1/delta_phi/sqrt(3));
end


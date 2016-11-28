function [ S_f ] = NaturalImages(cs, sv )
% Natural images 

fun_S = @(fr,ft) av(pi*ft/2./fr,sv)./fr.^3
integrand_fun = @(fr,ft) 2*pi*fr.*fun_S(fr,ft)
constant = integral2(integrand_fun,0,Inf,3/2/pi,240/2/pi)

S_f = @(fr,ft) cs*cs/constant*fun_S(fr,ft)


end


function [ S_f ] = FlatSpectrum(cs)
% Natural images 
max_ft = 300.0
fun_S = @(fr,ft) 1
integrand_fun = @(fr,ft) 2*pi*fr.*fun_S(fr,ft)
constant = integral2(integrand_fun,0,max_ft,3/2/pi,240/2/pi)

S_f = @(fr,ft) cs*cs/constant*fun_S(fr,ft)


end


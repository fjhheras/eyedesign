function information_rate = InformationRate( D, Ntu, S_f,f,Deltaphi, alias_fraction)
% Information of an eye with given characteristics per Hz and per sr

min_ft = 1.0; %Lower limit of time frequencies - Hz
max_ft = 300.0; %Upper limit of time frequencies - Hz
min_fr = 0;1/2/pi; %Lower limit of spatial frequencies, one cycle per 2pi - rad-1

d = 1.9; %rhabdomere diameter - um
lambda = 0.5; %wavelenght of yellow light - um
ImpulseWidth = 0.001; %0.001 sec -> for a fall of one order of magnitude at 170Hz 0.01 sec -> bslow 0.0001 sec -> bfast
LatencyWidth = 0.0014; %0.0014 sec -> for a corner frequency of 71Hz 0.0025 -> scan(2) 0.01 -> slow 0.02 -> slow2 0.04 -> slow3 0.1 -> superslow 0.00001 -> nolimit
MicrovilliRecycleTime = 0.02; %0.02; %sec

Deltarho2 =  lambda*lambda/D/D + d*d/f/f; %Normal one
%Deltarho2 =  1.5*lambda*lambda/D/D
coeffpoisson = Ntu/2/MicrovilliRecycleTime;

fb = 1/Deltaphi/sqrt(3); %array sampling frequency (rad-1)
intsampled = pi/2*fb*fb; %half a circle area of radius fb
%% Filters
lens_filter = @(fr) exp(-2*pi*pi/log(2)/4*(fr.*fr)*Deltarho2);
signal_filter = @(ft) 1./power(1+(2*pi*ImpulseWidth*ft).*(2*pi*ImpulseWidth*ft),3.12)*1./power(1+(2*pi*LatencyWidth*ft).*(2*pi*LatencyWidth*ft),2);
%% Signal
signal_power =  @(fr,ft) S_f(fr,ft).*lens_filter(fr).*signal_filter(ft);
signal_integrand = @(fr,ft) 2*pi*fr.*signal_power(fr,ft);
signal_variance = integral2(signal_integrand,min_fr,fb,min_ft,max_ft,'RelTol', 5e-4, 'AbsTol',1);
%% Noise
% Photon noise power. It depends on 2/coeffpoisson because we are in
% contrast. The extra factor two is because at the top of saturation SNR is
% half of what a poisson process with the same parameter
photon_noise_power = @(ft) 1./power(1+(2*pi*ImpulseWidth*ft).*(2*pi*ImpulseWidth*ft),3.12)/intsampled*2/coeffpoisson;
photon_noise_integrand = @(ft) 2./power(1+(2*pi*ImpulseWidth*ft).*(2*pi*ImpulseWidth*ft),3.12)*2/coeffpoisson;
photon_noise_variance = integral(photon_noise_integrand,min_ft,max_ft,'RelTol', 5e-4, 'AbsTol',1);

% Internal noise has the power spectrum shape of the signal, because we assume
% whitening before noise. It has the same power of photon noise.
internal_noise = @(fr,ft) 0; 
%internal_noise = @(fr,ft) signal_power(fr,ft)*photon_noise_variance/signal_variance; 

% A fraction (alias_fraction) of aliasing is considered noise 
aliased_power = @(fr,ft) alias_fraction.*signal_power(2*fb-fr,ft);

%% Information rates
H_integrand = @(fr,ft) fr.*log2( 1 + signal_power(fr,ft) ./ (photon_noise_power(ft) + internal_noise(fr,ft) + aliased_power(fr,ft))  );
information_rate = 2*sqrt(3)*integral2(H_integrand,min_fr,fb,min_ft,max_ft,'RelTol', 1e-3, 'AbsTol',1000);

end


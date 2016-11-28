D=30
Ntu = 100000
S_f=Sf
f = 60
Deltaphi = 0.01
alias_fraction=1

ft=logspace(0,3)

min_ft = 1; %Hz
max_ft = 300; %Hz
min_fr = 1/2/pi; %rad-1

d = 1.9; %rhabdomere diameter - um
lambda = 0.5; %wavelenght of yellow light - um
ImpulseWidth = 0.001; %0.001 sec -> for a fall of one order of magnitude at 170Hz 0.01 sec -> bslow 0.0001 sec -> bfast
LatencyWidth = 0.0014; %0.0014 sec -> for a corner frequency of 71Hz 0.0025 -> scan(2) 0.01 -> slow 0.02 -> slow2 0.04 -> slow3 0.1 -> superslow 0.00001 -> nolimit
MicrovilliRecycleTime = 0.02; %0.02; %sec
Deltarho2 =  lambda*lambda/D/D + d*d/f/f; %Normal one
cutoff_lens = 1/sqrt(Deltarho2); %kinda cutoff given by the lens. Used in aliasing
Deltatimp = 0.32/ImpulseWidth;
coeffpoisson = Ntu/2/MicrovilliRecycleTime;

fb = 1/Deltaphi/sqrt(3); %array sampling frequency (rad-1)
intsampled = pi/2*fb*fb; %half a circle area of radius fb

lens_filter = @(fr) exp(-2*pi*pi/log(2)/4*(fr.*fr)*Deltarho2);
signal_filter = @(ft) 1./power(1+(2*pi*ImpulseWidth*ft).*(2*pi*ImpulseWidth*ft),3.12)*1./power(1+(2*pi*LatencyWidth*ft).*(2*pi*LatencyWidth*ft),2);
signal_power =  @(fr,ft) S_f(fr,ft).*lens_filter(ft).*signal_filter(fr);
signal_integrand = @(fr,ft) 2*pi*fr.*signal_power(fr,ft);
signal_variance = integral2(signal_integrand,min_fr,fb,min_ft,max_ft);

photon_noise_power = @(ft) 1./power(1+(2*pi*ImpulseWidth*ft).*(2*pi*ImpulseWidth*ft),3.12)/intsampled*2/coeffpoisson;
photon_noise_integrand = @(ft) 2./power(1+(2*pi*ImpulseWidth*ft).*(2*pi*ImpulseWidth*ft),3.12)*2/coeffpoisson;
photon_noise_variance = integral(photon_noise_integrand,min_ft,max_ft);
photon_noise_variance/signal_variance
internal_noise = @(fr,ft) signal_power(fr,ft)*photon_noise_variance/signal_variance %internal noise power matching photon noise

aliased_power = @(fr,ft) alias_fraction.*signal_power(2*fb-fr,ft)

H_integrand = @(fr,ft) 2*pi*fr.*log( 1 + signal_power(fr,ft) ./ (photon_noise_power(ft) + internal_noise(fr,ft) + aliased_power(fr,ft))  )

information_rate = integral2(H_integrand,min_fr,fb,min_ft,max_ft)


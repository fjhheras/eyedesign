function H = InfG2( D, Ntu, N, Nt, fx,ft,S,f,Deltaphi,alias)

% Information of an eye with given characteristics
% Internal noise is supposed proportional to the shot noise so Nn not used
%


d = 1.9; %um
lambda = 0.5; %wave-lenght of yellow light in um

ImpulseWidth = 0.001; %0.001 sec -> for a fall of one order of magnitude at 170Hz 0.01 sec -> bslow 0.0001 sec -> bfast
LatencyWidth = 0.0014; %0.0014 sec -> for a corner frequency of 71Hz 0.0025 -> scan(2) 0.01 -> slow 0.02 -> slow2 0.04 -> slow3 0.1 -> superslow 0.00001 -> nolimit
ResolutionTime = 0.02; %0.02; %sec

Deltarho2 =  lambda*lambda/D/D + d*d/f/f; %Normal one
%Deltarho2 =  1.5*1.5*lambda*lambda/D/D; %Mod Stavenga

%%%% Temporal

%if (alias==1)
%    Deltarho2 =  lambda*lambda/D/D;
%else
%    Deltarho2 =  5*(lambda*lambda/D/D + d*d/f/f);
%end

%%%%



Deltatimp = 0.32/ImpulseWidth;
coeffpoisson = Ntu/2/ResolutionTime;

%fa = 1/Deltaphi/2;
fb = 1/Deltaphi/sqrt(3); %sampling frequency
%ftt = (power(10,1/3.12)-1)/2/pi/ImpulseWidth;
fD = 1/sqrt(Deltarho2);

if min(fb,fD) > fx(N)
    fb
    fx(N)
    pause;
end


nS = zeros(N,Nt);
nnS = zeros(N,Nt);
Snoise=zeros(N,Nt);

intsampled = 2*2*pi*pi*fb*fb*pi/2; %half a circle area of radius 2*pi*fb

for ir=1:N
    if (fx(ir) < min(fD,fb)) %it was max ?
        sum = 0;
        for it = 1:Nt
            nS(ir,it) = S(ir,it)*exp(-2*pi*pi/log(2)/4*(fx(ir)*fx(ir))*Deltarho2); %lens filtering
            nnS(ir,it) = nS(ir,it)*1/power(1+(2*pi*ImpulseWidth*ft(it))*(2*pi*ImpulseWidth*ft(it)),3.12)*1/power(1+(2*pi*LatencyWidth*ft(it))*(2*pi*LatencyWidth*ft(it)),2); %photoreceptor filtering ,24 -> scan
            Snoise(ir,it) = 1/power(1+(2*pi*ImpulseWidth*ft(it))*(2*pi*ImpulseWidth*ft(it)),3.12)/intsampled;
        end

    end
end

Ha=[0];
nnStotal=0;


sum2=0;
sumS=0;
sumN=0;

for ir = 2:N
    if (fx(ir) < min(fb,fD)) %it was fb
       sum = nnS(ir,1)/2;
       ssum=S(ir,1)/2;
       sumn=Snoise(ir,1)/2;
       for it = 2:Nt
          sum = sum + nnS(ir,it);
          ssum = ssum + S(ir,it);
          sumn = sumn + Snoise(ir,it);
       end
       
       sum2 = sum2 + (2*pi)*pi*fx(ir)*sum; %power of filtered signal
       sumS = sumS + (2*pi)*pi*fx(ir)*ssum;  %power of natural signal (not used)
       sumN = sumN + (2*pi)*pi*fx(ir)*sumn; %power of shot noise
end
end

%sum2=sum2/sumS*0.4*0.4;
%nnS=nnS./sumS*0.4*0.4; %normalizing everything to the same contrast as sampled

%noise = sum2*Nn/intsampled/ftt;
%

noise = 0; %internal noise proportional to noise in this case
Snoise = Snoise + nnS/sum2*sumN; %adding internal noise


%%% aliasing
%noise2=zeros(Nt);

if (alias==1 && fD > fb)

    sum3=zeros(Nt);
    for it = 1:Nt
        sum3(it)=(ceil(fb*2*pi)-fb*2*pi)*(2*pi)*pi*fx(ceil(fb*2*pi)-1)*nnS(ceil(fb*2*pi)-1,it);
        for ir = ceil(fb*2*pi):ceil(fD*2*pi)
         sum3(it) = sum3(it) + (2*pi)*pi*fx(ir)*nnS(ir,it);
        end
    
        if it==1
            sum3(it)=sum3(it)/2;
        end
    end
    noise2 = (sum3)./intsampled;

else
    noise2=zeros(Nt);        
end

%%% end-aliasing
%noise2=zeros(Nt);

%y=(1+2*real(1./(exp(-i*2*pi*ResolutionTime*ft).*(1-i*2*pi*ResolutionTime*ft)-1))); %correction

Tnoise=0;
Tsignal=0;
for ir = 2:N
    if (fx(ir) < min(fb,fD))
       sum = log(1+(nnS(ir,1)*coeffpoisson/4)/(Snoise(ir,1)/2+(noise+noise2(1))*coeffpoisson/4))/2;
       for it = 2:Nt
          sum = sum + log(1+(nnS(ir,it)*coeffpoisson/4)/(Snoise(ir,it)/2+(noise+noise2(it))*coeffpoisson/4));
          %(nnS(ir,it)*coeffpoisson/nSnorm)/(Snoise(ir,it))
          %Tnoise=Tnoise+(Snoise(ir,it)+(noise+noise2(it))*coeffpoisson);
          %Tsignal=Tsignal+nnS(ir,it)*coeffpoisson;
       end
       
       %if (fx(ir) > fa)
       %    sum = sum*(1-6*acos(1/fx(ir)/2/Deltaphi)/pi); 
       %end
       
       Ha(ir) = pi*(2*pi*fx(ir))*sum;
end
end
%H = H * Deltax*Deltay*Deltat;

%StoN=sqrt(Tsignal/Tnoise)

H = trapz(Ha);

end


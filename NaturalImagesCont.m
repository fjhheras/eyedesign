function [ S ] = NaturalImagesCont(  N, Nt,cs, sv )
% Natural images 

fun_S = @(fr,ft) av(pi*ft/2./fr,sv)./fr.^3
integrand_fun = @(fr,ft) 2*pi*fr.*fun_S(fr,ft)
constant = integral2(integrand_fun,0,Inf,3/2/pi,240/2/pi)


W=2*pi %rads
%Wy=2*pi
Wt=1 %secs
%Wr=1/sqrt(1/W/W+1/W/W)

Deltax = 1/W;
%Deltay = 1/Wy;
Deltat = 1/Wt;

nx=0:N-1; % l -|m|
fx = nx./W;

nt=0:Nt-1; %-Nt/2+1:Nt/2
ft = nt./Wt;

ix =1
for it = 1:Nt
        S(ix,it)= 0;
end  
    
for ix=2:N
    for it = 1:Nt
        S(ix,it)= cs*cs/constant*fun_S(fx(ix),ft(it))*Deltax*Deltax*Deltat;
    end    
end

end


function [ S, fx, ft, Deltax, Deltat ] = NatImagbigSHt( N, Nt, cs, sv )
% Natural images with spherical harmonics
% N = number of spatial sampling points, determines max spatial frequency
% Nt = number of temporal sampling points, determines max temporal frequency
% cs = contrast 
% sv = angular velocity in rad/s
% Returns: double S[N][Nt]
% Where S[ix][it]/Deltax/Deltay/Deltat is the power spectrum of natural images

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
ftx = (-500000:500000)/Wt;

c1=0;

for ix = 1:240
    for iy = 1:240
            ir2 = (ix-1)*(ix-1)+(iy-1)*(iy-1);
            if (ir2>9 && ir2 < 57600) 
                fr2 = (fx(ix)*fx(ix)+fx(iy)*fx(iy));
                c1= c1 + 1/fr2;
                if (ix>1)
                    c1= c1 + 1/fr2; % m negative values also count
                end
            end
    end
end
    
cs

iy=1;

below_threshold = 1;
for ix=2:N
    ix
    fr = fx(ix);
    
    %At each step we calculate cv, a variable that must be chosen so sum av(v)dv = 0 (see appendix of van Hateren, 1992)
	% We calculate it using the sum given by eq 9, but we extend it to more values of ft.
	% We save some computing time by realising that, at high spatial frequencies, when the program runs out of
	% temporal frequencies, estimated cv grows above sv/2, while it should asymptote.
    
    if below_threshold      
        cv = 1/( pi/2*sum(av(pi*ftx/2./fr,sv)./Wt./fr))
        if cv > sv/2
            cv = sv/2
            below_threshold = 0
        end
    end
    

    for it = 1:Nt
        S(ix,it)= pi*cs*cs*cv*av(pi*ft(it)/2/fr,sv)/(2*c1*(1+cs*cs)*fr*fr*fr);
    end    
    
end

S(1,1)=Wt/(1+cs*cs);
for it=2:Nt
   S(1,it)=0; 
end

S = S/S(1,1);

end


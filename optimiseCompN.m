function [Da, Na, Ra, Ha, volume_a] = optimiseCompN(c)

%clear all;
Sf = NaturalImages(0.4,1)
%Sf = FlatSpectrum(0.4)
LoadConstants
LoadCalliphoraV0

DPoints=100;
NPoints=200;
invpoints=20

%Darray = linspace(15,100,DPoints);
Darray = logspace(1,2,DPoints);
cagada=0;
%Narray2 = linspace(3,6,50)
inva = logspace(-1,2.1,invpoints) * V0
%Harray2 = zeros(20,50)


Da=zeros(size(inva));
Na=zeros(size(inva));
Ra=zeros(size(inva));
Ha=zeros(size(inva));
Hmaxa = zeros(size(inva));
volume_a = zeros(size(inva));
volume_debugging_a = zeros(size(inva));


%c=0.7e-10
%c=1e-10
%c=10e-10
%c=0.0
%c = c*V0


invi=0;
for inv = inva
invi=invi+1;
inv

%Narray = logspace(2.5,6.0,NPoints);
maxpossibleL = nthroot(inv / (4*pi/3),3);
Narray = logspace(2.5,log10(maxpossibleL/k),NPoints);
Hmax=0;
Dmax =-1;
Nmax =-1;
for iD=1:DPoints
    for iN = 1:NPoints
        H(iD,iN) = InfCompN(Darray(iD),Narray(iN), Sf, inv, c);
        if H(iD,iN)>Hmax
            Hmax = H(iD,iN);
            Dmax=Darray(iD);
            Nmax=Narray(iN);
        end
    end
end

options = optimset('TolX',.001);
fH = @(x) -InfCompN( x(1), x(2),Sf, inv,c );
[y,H] = fminsearch(fH,[Dmax, Nmax], options);
Ha(invi) = -H;
Da(invi) = y(1);
Na(invi) = y(2);

if (InfCompN( Dmax, Nmax, Sf, inv,c ) > InfCompN( y(1), y(2),Sf, inv, c))
    cagada = cagada+1
    y(1)=Dmax;
    y(2)=Nmax;
end

R = K2R(inv/(4*pi/3),y(1),y(2),c,f_number,k,6, n3);
volume_a(invi) = inv - V_equiv_energy_cost(c,6,y(1),y(2),R ); %6*c*4*pi*y(2)/y(1)/y(1)/(2*sqrt(3)/3)*R*R
volume_debugging_a(invi) = R2K(R,y(1),y(2),0,f_number,k,6,n3)*4*pi/3; %(R*R*R-(R-L)*(R-L)*(R-L))/K0

%% DEBUGGING - Uncomment this to test K2R with RK2
should_be_one = K2R(R2K(R,y(1),y(2),0,f_number,k,6,n3),y(1),y(2),0,f_number,k,6,n3)/R
%%

Ra(invi) = R;
Hmaxa(invi) = Hmax; 
end


pa=Da.*Da./Ra;
La = Na*k;% + cone_length(Da,f_number,n3);



function [Da, Na, Ra, Ha, volume_a] = optimiseCompN(c, inva)

%clear all;
Sf = NaturalImages(0.4,1)
%Sf = FlatSpectrum(0.4)
LoadConstants
%LoadCalliphoraV0

%% Parameters for fast test
DPoints=40;
NPoints=15;

%% Parameter values for long computation
%DPoints=200;
%NPoints=50;

%% Building arrays
Darray = logspace(0.5,2.5,DPoints);


Da=zeros(size(inva));
Na=zeros(size(inva));
Ra=zeros(size(inva));
Ha=zeros(size(inva));
Hmaxa = zeros(size(inva));
volume_a = zeros(size(inva));
volume_debugging_a = zeros(size(inva));

%% Loop
invi=0;
errors=0;
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
    
    % Taking the best among the ones tried, we use it
    % as starting point of the optimization
    options = optimset('TolX',.001);
    fH = @(x) -InfCompN( x(1), x(2),Sf, inv,c );
    [y,H] = fminsearch(fH,[Dmax, Nmax], options);
    Ha(invi) = -H;
    Da(invi) = y(1);
    Na(invi) = y(2);

    if (InfCompN( Dmax, Nmax, Sf, inv,c ) > InfCompN( y(1), y(2),Sf, inv, c))
        errors = errors+1
        y(1)=Dmax;
        y(2)=Nmax;
    end

R = K2R(inv/(4*pi/3),y(1),y(2),c,f_number,k,6, n3);
volume_a(invi) = inv - V_equiv_energy_cost(c,6,y(1),y(2),R ); 
volume_debugging_a(invi) = R2K(R,y(1),y(2),0,f_number,k,6,n3)*4*pi/3; 

%% DEBUGGING - Uncomment this to test K2R with RK2
should_be_one = K2R(R2K(R,y(1),y(2),0,f_number,k,6,n3),y(1),y(2),0,f_number,k,6,n3)/R
%%

Ra(invi) = R;
Hmaxa(invi) = Hmax; 
end


%pa=Da.*Da./Ra;
%La = Na*k;% + cone_length(Da,f_number,n3);
%pause

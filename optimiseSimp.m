function [Da, Na, Ra, Ha, volume_a] = optimiseSimp(c, inva)

Sf = NaturalImages(0.4,1)
LoadConstants

%% Parameters for fast test
%NPoints=15;

%% Parameter values for long computation
NPoints=80;

%% Building arrays
Da=zeros(size(inva));
Na=zeros(size(inva));
Ra=zeros(size(inva));
Ha=zeros(size(inva));
Hmaxa = zeros(size(inva));
volume_a = zeros(size(inva));
volume_debugging_a = zeros(size(inva));

%% Loop
invi=0;

for inv = inva
    invi=invi+1;
    inv

    % Calculating an array of num microvilli to sample
    maxpossibleL = nthroot(inv / (4*pi/3),3);
    Narray = logspace(2.5,log10(maxpossibleL/k),NPoints);


    % Sampling a matrix of values to find a starting point for optimization
    Hmax=0;
    Dmax = -1;
    Nmax = -1;
    for iN = 1:NPoints
        H(iN) = InfSimp(Narray(iN), Sf, inv, c);
        if H(iN)>Hmax
            Hmax = H(iN);
            Nmax=Narray(iN);
        end
    end
    
    % Taking the best among the ones tried, we use it
    % as starting point of the optimization
    options = optimset('TolX',1e-7);
    fH = @(x) -InfSimp(x, Sf, inv, c);
    [y,H] = fminsearch(fH, Nmax, options);
    Ha(invi) = -H;
    Da(invi) = y;

R = K2R(inv/(4*pi/3), y(1), y(2), c, f_number, k, b, n3);


volume_a(invi) = inv - V_equiv_energy_cost(c,b,y(1),y(2),R ); 
volume_debugging_a(invi) = R2K(R,y(1),y(2),0,f_number,k,b,n3)*4*pi/3; 

%% DEBUGGING - Uncomment this to test K2R with RK2
%should_be_one = K2R(R2K(R,y(1),y(2),0,f_number,k,b,n3),y(1),y(2),0,f_number,k,b,n3)/R
%%

Ra(invi) = R;
Hmaxa(invi) = Hmax; 
end

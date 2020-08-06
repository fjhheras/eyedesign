function [Da, Na, Ra, Ha, volume_a] = optimiseSimp(c, inva)

% 1 or 3 for comparison against comp
b = 1;

Sf = NaturalImages(0.4,1)
LoadConstants

%% Parameters for fast test
%NPoints=11;

%% Parameter values for long computation
NPoints=25;

%% Building arrays
Da=zeros(size(inva));
Na=zeros(size(inva));
Ra=zeros(size(inva));
Ha=zeros(size(inva));
Hmaxa = zeros(size(inva));
volume_a = zeros(size(inva));

%% Loop
invi=0;
for inv = inva
    invi=invi+1
    
    % Calculating an array of num microvilli to sample
    maxpossibleL = nthroot(inv / (4*pi/3),3);
    Narray = logspace(1.5,log10(maxpossibleL/k),NPoints);


    % Sampling a matrix of values to find a starting point for optimization
    Hmax = 0;
    Dmax = -1;
    Nmax = -1;
    for iN = 1:NPoints
        H(iN) = InfSimp(Narray(iN), Sf, inv, c, b);
        if H(iN)>Hmax
            Hmax = H(iN);
            Nmax = Narray(iN);
        end
    end
    
    % Taking the best among the ones tried, we use it
    % as starting point of the optimization
    options = optimset('TolX',1e-7);
    fH = @(x) -InfSimp(x, Sf, inv, c, b);
    [y,H] = fminsearch(fH, Nmax, options);
    Ha(invi) = -H;
    Na(invi) = y;

    L = k * Na(invi);
    R = R_in_simple(inv, L, d, c, Na(invi), b);
    f = R/n3; % focal length
    Da(invi) = f/f_number; % Lens diameter from f-number and focal length
    volume_a(invi) = (4*pi/3) * (R + L)^3; 
    
    
    %% DEBUGGING - Uncomment this to test
    volume_a_check = inv - V_equiv_energy_cost(c, b, d, Na(invi), R );
    should_be_one = volume_a(invi)/volume_a_check
    %%
    
    Ra(invi) = R;
    Hmaxa(invi) = Hmax; 
end

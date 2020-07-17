%clear all;
Sf = NaturalImages(0.4,1)
LoadConstants
LoadCalliphoraV0

DPoints=20;
NPoints=80;


Darray = linspace(12,45,DPoints);
Narray = logspace(2,5.2,NPoints);

inv = V0

%c=0.7e-10
%c=1e-10
%c=10e-10
c=0.0

H=[];
for iD=1:DPoints
    for iN = 1:NPoints
        H(iD,iN) = InfCompN(Darray(iD), Narray(iN), Sf, inv, c);
    end
end

%% Log plot in number of microvilli
surf(log10(Narray), Darray,H)

%% Linear plot in both number of microvilli and lens diameter
%surf(Narray, Darray,H)

pause


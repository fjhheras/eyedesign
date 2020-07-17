%clear all;
Sf = NaturalImages(0.4,1)
LoadConstants
LoadCalliphoraV0

DPoints=20;
NPoints=80;


H=[];
Darray = linspace(12,45,DPoints);
Narray = logspace(2,5.2,NPoints);

inv = V0

%c=0.7e-10
%c=1e-10
%c=10e-10
c=0.0

for iD=1:DPoints
    for iN = 1:NPoints
        H(iD,iN) = InfCompN(Darray(iD),Narray(iN), Sf, inv, c);
    end
end


surf(log10(Narray), Darray,H)
%surf(Narray, Darray,H) %if we fix cost by fixing volume, dividing by cost again is not necessary (and it's wrong cos R is no longer ctant)

pause


%clear all;
Sf = NaturalImages(0.4,1)
LoadConstants
LoadCalliphoraV0

DPoints=40;
NPoints=150;

Darray = linspace(12,45,DPoints);
Narray = logspace(2,5.1,NPoints);

%c=0.7e-10
%c=1e-10
%c_a=[0.0 0.05 0.1 0.5]
c=0.13

% 4e9 microvilli in the whole Calliphora eye
inv = V0 + c * 4e9

H=[];
for iD=1:DPoints
    for iN = 1:NPoints
        H(iD,iN) = InfComp(Darray(iD), Narray(iN), Sf, inv, c, 'fly');
    end
end
filename = sprintf('%s_%0.3e','info3d',c)
csvwrite(filename, H)

%% Log plot in number of microvilli
surf(log10(Narray), Darray,H)
xlabel('log number microvilli')
ylabel('Lens diameter (um)')
zlabel('Information rate (bit Hz/sr)')
%% Linear plot in both number of microvilli and lens diameter
%surf(Narray, Darray,H)
    


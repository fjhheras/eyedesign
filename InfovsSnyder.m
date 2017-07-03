Da = linspace(10,50)
R = 1200
[H, Ha, Hs] = InformationRates(Da,R)
plot(Da,H,'r')
hold on
plot(Da,Ha,'r--')
plot(Da,Hs,'k')

R = 400
[H, Ha, Hs] = InformationRates(Da,R)
plot(Da,H,'b')
hold on
plot(Da,Ha,'b--')
plot(Da,Hs,'k')

xlabel('Lens diameter (um)')
ylabel('Information rate bits/s/sr')

function [H_ST, H_ST_with_aliasing, H_HowardSnyder] = InformationRates(Da, R)

    cs = 0.4
    sv = 1.0 %0.3
    Sf = NaturalImages(cs, sv)

    MicrovilliRecycleTime = 0.02; %0.02; %sec
    NumberMicrovilli = 60000;
    F_number = 2
    N = NumberMicrovilli / MicrovilliRecycleTime


    H_ST=[]
    H_ST_with_aliasing=[]
    H_HowardSnyder=[]
    for D = Da
        H_ST=[H_ST,InformationRate(D,NumberMicrovilli,Sf,D*F_number,D/R,0)];
        H_ST_with_aliasing=[H_ST_with_aliasing,InformationRate(D,NumberMicrovilli,Sf,D*F_number,D/R,1)];
        H_HowardSnyder = [H_HowardSnyder, snyderH( R, D, F_number, cs , N )]
    end
end









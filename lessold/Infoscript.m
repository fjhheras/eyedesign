cs = 0.4
sv = 1.0 %0.3
Sf = NaturalImages(cs, sv)
Da = linspace(10,60)

R=1200

H=[]
for D = Da
    H=[H,InformationRate(D,500000,Sf,D*2,D/R,0)];
end

Ha=[]
for D = Da
    Ha=[Ha,InformationRate(D,500000,Sf,D*2,D/R,1)];
end

plot(Da,[H;Ha],'--') 

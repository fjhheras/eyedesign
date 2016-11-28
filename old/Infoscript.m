N =3000
Nt =250
cs = 0.4
sv = 0.3
[ S,fx,ft ] = NaturalImagesCont(  N, Nt,cs, sv )

Da = linspace(1,60)
H=[]
R=1200

for D = Da
    H=[H,InfG2(D, 500000, N, Nt, fx,ft,S,D*2,D/R,0)];
end

Ha=[]
for D = Da
    Ha=[Ha,InfG2(D, 500000, N, Nt, fx,ft,S,D*2,D/R,1)];
end

plot(Da,[H;Ha])
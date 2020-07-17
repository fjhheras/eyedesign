LoadCalliphoraV0
%c = [0 1e-12 3e-12 1e-11 0.3e-10 1e-10 3e-10 1e-9]*V0
c = [0 0.04 0.08 0.16 0.32 0.64]
numc = 6
Naa = zeros(numc, 20)
Daa = zeros(numc, 20)
Haa = zeros(numc, 20)
volume_aa = zeros(numc, 20)
parfor i = 1:numc
    [Da, Na, Ra, Ha, volume_a] = optimiseCompN(c(i))
    Naa(i,:) = Na
    Daa(i,:) = Da
    Haa(i,:) = Ha
    Raa(i,:) = Ra
    volume_aa(i,:) = volume_a
end

um3_to_mm3_per_sr = 1e-9 / (4 * pi)

figure(1)

LoadConstants
for i = 1:numc
    loglog(volume_aa(i,:)*um3_to_mm3_per_sr, Naa(i,:)*k)
    hold on
end
xlabel('Volume mm3/sr')
ylabel('Rh length (um)')

figure(2)

for i = 1:numc
    loglog(volume_aa(i,:)*um3_to_mm3_per_sr, Haa(i,:))
    hold on
end
xlabel('Volume mm3/sr')
ylabel('Information rate')


figure(3)

for i = 1:numc
    loglog(volume_aa(i,:)*um3_to_mm3_per_sr, Daa(i,:).* Daa(i,:) ./ Raa(i,:) )
    hold on
end
xlabel('Volume mm3/sr')
ylabel('Eye parameter (um)')
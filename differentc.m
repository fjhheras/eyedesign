% Extra cost of photoreceptor material in equivalent volume
% Units um3 / microvilli
% Note: in paper estimated 6e-11 - 3.5e-10 in mm3/microvilli 
% which is 0.06 - 0.35
c = [0 0.04 0.08 0.16 0.32 0.64]
numc = length(c)

% Goal: 3e‑6 to 1.5e2 mm3/sr, 3.8e4 to 1.9e12 um3
% Avoiding low range and compensating for cost
inva = logspace(6.5, 12.5, 20)
num_inv = length(inva)

Naa = zeros(numc, num_inv);
Daa = zeros(numc, num_inv);
Haa = zeros(numc, num_inv);
volume_aa = zeros(numc, num_inv);
%parfor
parfor i = 1:numc
    [Da, Na, Ra, Ha, volume_a] = optimiseCompN(c(i), inva)
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
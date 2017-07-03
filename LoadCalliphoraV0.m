%%Estimating Calliphora eye volume

Deltaphi0 = 2*pi/180;
R0 = 31/(Deltaphi0);
L0 = 250 + 60; %um
K0 = 3*L0*R0*R0 - 3*L0*L0*R0+L0*L0*L0;
V0 = K0*4*pi/3
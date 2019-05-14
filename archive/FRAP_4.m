% A=0.98;
% C=0.7;
% Tau=2.27e3;
A = 0.7;
C = 0.15;
Tau = 1.5e4;

icalc = A.*exp(-Tau./(2.*tt)).*(besseli(0,Tau./(2.*tt))+besseli(1,Tau./(2.*tt)))+C;
plot(tt, ii)
hold on
plot(tt,ii1)
plot(tt, icalc)
hold off
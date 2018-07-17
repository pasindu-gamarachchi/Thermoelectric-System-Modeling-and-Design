function Qh = TEBlock(Th,Tc)

 global cvl hxw

% cvl = cvl; %[m]
k = 1.225631; %2.2638; % value for Nicks Model = ; 

hl = 4.90*10^-3; % [m]


A = hxw*cvl;
Rth = hl/(k*A);

Qh = (Th -Tc)/(Rth);
%% Power Generation in a TE Couple, with Power output of 
% two legs seperated
 function [Qh P]  = Unicouple(Th, Tc)


% Naming Convention
% property/leg/otherconsideration   no slashes

%  Th = 491.48;  % 300 400 500 600 ];% [C] Hot side temperature
%m = 94;
% RL = 0:0.001:0.8;
%  for i= 1:5;
% Th= 500;
% Tc = 94;   % [C] Cold side temperature
hn = 0.002 ; % [m] n-leg height 
ln = 0.0018; %  [m] n-leg length
tn = 0.0018; %  [m] n-leg thickness
hp = 0.002; %  [m] p-leg height
lp = 0.0018; % [m] p-leg length
tp =0.0018 ; % [m] p-leg thickness

%% Temperature dependent properties/ Integrals of functions with respesct to T
%  Resistivity N
syms x
pnf = int(-1.16561436E-14*x^3 + 5.03594346E-12*x^2 + 6.47260002E-09*x + 6.50670403E-06);

%% Seebeck N
syms x
snf = int(5.48885661E-14*x^3 + 2.22674400E-10*x^2 - 2.66412731E-07*x - 1.33513735E-04);

%% Conductivity N
syms x
knf = int(1.18722952E-08*x^3 - 4.65123739E-06*x^2 - 7.31176479E-04*x + 3.53715861E+00);

%% Resisitivity P
syms x
ppf = int(-1.26642460E-15*x^3 - 6.37965919E-12*x^2 + 1.55655066E-08*x + 1.03160255E-05);

%% Seebeck P
 syms x
 spf= int( 1.50803298E-13*x^3 - 3.37285775E-10*x^2 + 3.01002672E-07*x + 1.23997961E-04);

 %% Conductivity P
 
 syms x
 kpf = int( -1.79754155E-16*x^6 + 3.96173694E-13*x^5 - 3.42301654E-10*x^4 + 1.48422824E-07*x^3 - 3.25292137E-05*x^2 + 2.66399769E-03*x + 3.06261874E+00);
 
 %%
syms kn kp sn sp pp pn x % k: Conductivity, s: Seebeck, p: resistivity, p: p-leg, n:n-leg


kn = symfun( knf, x); % Temperature dependent k for n-leg
kp = symfun (kpf, x); % Temperature dependent k for p-leg
sn = symfun (snf, x); % Temperature dependent s for n-leg
sp = symfun (spf, x); % Temperature dependent k for p-leg
pp = symfun (ppf, x); % Temperature dependent k for p-leg
pn = symfun (pnf, x); % Temperature dependent k for p-leg
%% Integral Averages 

knm = (kn(Th) - kn (Tc))./(Th-Tc);
knm = double(knm);
kpm = (kp(Th) - kp (Tc))./(Th-Tc);
kpm = double(kpm);
snm =  (sn(Th) - sn(Tc))./(Th-Tc);
snm = double(snm);
spm =  (sp(Th) - sp(Tc))./(Th-Tc);
spm = double(spm);
ppm =  (pp(Th) - pp(Tc))./(Th-Tc);
ppm = double(ppm);
pnm =  (pn(Th) - pn(Tc))./(Th-Tc);
pnm = double(pnm);


%% Calculations

rp = (ppm .* hp )./(lp.*tp);  % [Ohms] Resistance p-leg
rn = (pnm .* hn )./(ln.*tn); % [Ohms] Resistance n-leg
cp =  (kpm .* lp.*tp)./hp; %  Thermal Conductance of p-leg
cn = (knm .* ln.*tn)./hn; % Thermal Conductance of n-leg

CR =

Rt = rp +rn; % [Ohms] Resistance Total 
RL = Rt; % [Ohms] Load resistance is set to internal resistance
R = RL +Rt; % [Ohms] 
S =  spm -snm ; % Seebeck of couple
I = (S.*(Th -Tc))./(R);% [A] Current through the unicouple

%% P- Leg Power Output and efficiency

Qhp = spm.*Th.*I + cp.*(Th-Tc) - (0.5).*(I.^2).*rp;  % Heat input to hot side
Qcp = spm.*Tc.*I + cp.*(Th-Tc) + 0.5.*(I.^2).*rp; % Heat rejection from cold side
Pp = Qhp- Qcp ;% Power output per couple
Np = (Pp./Qhp).*100;


%% N- Leg Power Output and efficiency

Qhn = abs(snm).*Th.*I + cn.*(Th-Tc) - (0.5).*(I.^2).*rn ; % Heat input to hot side
Qcn = abs(snm).*Tc.*I + cn.*(Th-Tc) + 0.5.*(I.^2).*rn ;% Heat rejection from cold side
Pn = Qhn- Qcn; % Power output per couple
Nn = (Pn./Qhn).*100;

%% Total Power

P = Pp + Pn;
N = P./(Qhp +Qhn)*100;
%  end
Qh = Qhp + Qhn;
Pn;
Pp;
P;

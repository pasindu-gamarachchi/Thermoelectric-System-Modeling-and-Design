%% Power Generation in a TE Couple, with Power output of 
% two legs seperated
% Pasindu Gamarachchi
%
clc
clear 
close all
format long
tic

% Updated 06/09
Date = date
DBC_Dimensions = 'smaller dimensions from drawings';
TE_Properties = 'Half - Heusler + Bi2Te3'; % measured 24/05 by Nick
Leg_Dimensions = 'Dimensions of unicouple for DOE';
% Updated 06/15
% Alumina and Copper Properties Updated



% Naming Convention
% property/leg/otherconsideration   no slashes

%  Th = [200 300 400 500 600 ];% [C] Hot side temperature
elems =500; % m is the number of segments, Set m >1
numn = 2*elems +1;

Th = 600;
Tc = 20;   % [C] Cold side temperature

TLH = 2.4*(10^-3);
% Dimensions

% N - Leg
hn = (2.1056*10^-3) ;% (1.759*10^-3) ; % [m] n-leg height 
ln = 1.5*10^-3; %  [m] n-leg length
tn = 1.5*10^-3; %  [m] n-leg thickness
% P - Leg

hp = (2.1056*10^-3); %  [m] p-leg height
lp = 1.5*10^-3; % [m] p-leg length
tp = 1.5*10^-3 ; % [m] p-leg thickness

% N - Leg 2
hn2 = TLH - hn ; % [m] n-leg height 
ln2 = 1.5*10^-3; %  [m] n-leg length
tn2 = 1.5*10^-3; %  [m] n-leg thickness
% P - Leg 2

hp2 = TLH -hp ; %  [m] p-leg height
lp2 = 1.5*10^-3; % [m] p-leg length
tp2 = 1.5*10^-3 ; % [m] p-leg thickness




%  T3 - Copper
t3w = 1.93*10^-3;
t3l = 1.96*10^-3;
t3t = 0.2032*10^-3;

% T2 - Ceramic
t2w = 4.51*10^-3;
t2l = 2.26*10^-3;
t2t = 0.635*10^-3;
% it2t = t2t;

% T1 - Copper
t1w = 4.21*10^-3;
t1l = 1.96*10^-3;
t1t = 0.2032*10^-3;

% B3 - Copper
b3w = 8.50*10^-3;
b3l = 1.96*10^-3;
b3t = 0.2032*10^-3;

% B2- Ceramic
b2w = 8.81*10^-3;
b2l = 2.26*10^-3;
b2t = 0.635*10^-3;

% B1 - Copper
b1w = 4.07*10^-3;
b1l = 1.96*10^-3;
b1t = 0.2032*10^-3;

ws = warning('off', 'all');

C_pho = 10*(10^-6);

% Temperature Dependent properties

Notes = 'Material Properties Measured by Nick 05/24/2016'
syms T

% Material Properties Top Material
% Seebeck - N 
TNS = [ 20 ,50, 100, 200, 300, 400, 500, 600 ];
SN = [ -0.000135243, -0.000142471, -0.000153643, -0.000171607, -0.000187049, -0.000202244, -0.000211448, -0.000217053];
CNS = polyfit( TNS, SN, 6);
snf = int ( CNS(1)*T^6 + CNS(2)*T^5 + CNS(3)*T^4 + CNS(4)*T^3 + CNS(5)*T^2 + CNS(6)*T + CNS(7));

% Electrical - N 
TNP = [ 20 ,50, 100, 200, 300, 400, 500, 600 ];
PN = [ 6.56786E-06, 6.7325E-06, 7.03657E-06, 7.60784E-06, 8.0957E-06, 8.4766E-06, 8.6924E-06, 8.77723E-06];
CNP = polyfit( TNP, PN, 6);
pnf = int ( CNP(1)*T^6 + CNP(2)*T^5 + CNP(3)*T^4 + CNP(4)*T^3 + CNP(5)*T^2 + CNP(6)*T + CNP(7));

% Theramal - N 
TNK = [ 20 ,50, 100, 200, 300, 400, 500, 600 ];
KN = [ 5.308966362, 5.001685214, 4.609217214, 4.251736671, 4.069193415, 3.95510388, 3.962709849, 4.175676981];
CNK = polyfit( TNK, KN, 6);
knf = int ( CNK(1)*T^6 + CNK(2)*T^5 + CNK(3)*T^4 + CNK(4)*T^3 + CNK(5)*T^2 + CNK(6)*T + CNK(7));

% Seebeck - P 
TPS = [ 20 ,50, 100, 200, 300, 400, 500, 600 ];
SP = [7.39928E-05, 7.95446E-05, 8.87751E-05, 0.00010694, 0.000125184, 0.000142737, 0.000161509, 0.00017735 ];
CPS = polyfit( TPS, SP, 6);
spf = int ( CPS(1)*T^6 + CPS(2)*T^5 + CPS(3)*T^4 + CPS(4)*T^3 + CPS(5)*T^2 + CPS(6)*T + CPS(7));

% Electrical - P
TPP = [ 20 ,50, 100, 200, 300, 400, 500, 600 ];
PP = [1.97399E-06, 2.13811E-06, 2.46366E-06, 3.18329E-06, 4.01544E-06, 4.92E-06, 5.85437E-06, 6.7556E-06 ];
CPP = polyfit( TPP, PP, 6);
ppf = int ( CPP(1)*T^6 + CPP(2)*T^5 + CPP(3)*T^4 + CPP(4)*T^3 + CPP(5)*T^2 + CPP(6)*T + CPP(7));

% Thermal - P
TPK = [ 20 ,50, 100, 200, 300, 400, 500, 600 ];
KP = [7.533937719, 7.251196645, 6.711129425, 5.980450245, 5.438397484, 5.05121694, 4.773241165, 4.550860545];
CPK = polyfit( TPK, KP, 6);
kpf = int ( CPK(1)*T^6 + CPK(2)*T^5 + CPK(3)*T^4 + CPK(4)*T^3 + CPK(5)*T^2 + CPK(6)*T + CPK(7));

% Material 2
% Seebeck - N 
TNS2 = [ 26.4,	51.1,	74.9,	99.3,	124.7,	150.3,	174.9,	200.2,	223.9,	249.8 ];
SN2 = [ -1.90E-04,	-1.96E-04,	-2.04E-04,	-2.03E-04,	-2.05E-04,	-2.04E-04,	-2.00E-04,	-1.93E-04,	-1.86E-04,	-1.71E-04];
CNS2 = polyfit( TNS2, SN2, 6);
snf2 = int ( CNS2(1)*T^6 + CNS2(2)*T^5 + CNS2(3)*T^4 + CNS2(4)*T^3 + CNS2(5)*T^2 + CNS2(6)*T + CNS2(7));

% Electrical - N 
TNP2 = TNS2;
PN2 = [ 1.01E-05,	1.11E-05,	1.22E-05,	1.33E-05,	1.44E-05,	1.55E-05,	1.64E-05,	1.71E-05,	1.75E-05,	1.77E-05];
CNP2 = polyfit( TNP2, PN2, 6);
pnf2 = int ( CNP2(1)*T^6 + CNP2(2)*T^5 + CNP2(3)*T^4 + CNP2(4)*T^3 + CNP2(5)*T^2 + CNP2(6)*T + CNP2(7));

% Theramal - N 
TNK2 = TNS2;
KN2 = [1.26,	1.24,	1.23,	1.20,	1.22,	1.29,	1.30,	1.37,	1.47,	1.57];
CNK2 = polyfit( TNK2, KN2, 6);
knf2 = int ( CNK2(1)*T^6 + CNK2(2)*T^5 + CNK2(3)*T^4 + CNK2(4)*T^3 + CNK2(5)*T^2 + CNK2(6)*T + CNK2(7));

% Seebeck - P 
TPS2 = [ 25.7,	51.8,	75.5,	99.8,	125.4,	151,	175.5,	200.9,	224.4,	250.3];
SP2 = [0.000192,	0.000202,	0.000211,	0.000216,	0.000216,	0.000218,	0.000216,	0.000211,	0.000205,	0.000194 ];
CPS2 = polyfit( TPS2, SP2, 6);
spf2 = int ( CPS2(1)*T^6 + CPS2(2)*T^5 + CPS2(3)*T^4 + CPS2(4)*T^3 + CPS2(5)*T^2 + CPS2(6)*T + CPS2(7));

% Electrical - P
TPP2 = TPS2;
PP2 = [9.96E-06,	1.15E-05,	1.30E-05,	1.47E-05,	1.66E-05,	1.84E-05,	2.02E-05,	2.19E-05,	2.33E-05,	2.45E-05 ];
CPP2 = polyfit( TPP2, PP2, 6);
ppf2 = int ( CPP2(1)*T^6 + CPP2(2)*T^5 + CPP2(3)*T^4 + CPP2(4)*T^3 + CPP2(5)*T^2 + CPP2(6)*T + CPP2(7));

% Thermal - P
TPK2 = TPS2;
KP2 = [1.08,	1.05,	1.03,	1,	1,	1,	1.03,	1.07,	1.11,	1.18];
CPK2 = polyfit( TPK2, KP2, 6);
kpf2 = int ( CPK2(1)*T^6 + CPK2(2)*T^5 + CPK2(3)*T^4 + CPK2(4)*T^3 + CPK2(5)*T^2 + CPK2(6)*T + CPK2(7));





% Thomson N  - New code 
dsdtnf = symfun(diff(snf,T),T); % Derivative of seebeck with respect to T
thonf = symfun( (dsdtnf),T); % Temeperature dependent Thomson Coefficient
thnf =int(thonf); % Integral of temperature dependent Thomson Coefficient

% Thomson P  - New code 
dsdtpf = symfun(diff( spf,T),T); % Derivative of seebeck with respect to T
thopf = symfun((dsdtpf),T);   % Temeperature dependent Thomson Coefficient
thpf = int(thopf); % Integral of temperature dependent Thomson Coefficient
 
% Conductivity Copper 102 
TCK = [2.85, 27.85,	77.85, 127.85, 177.85, 227.85, 277.85 , 327.85, 377.85, 427.85,	477.85,	527.85,	577.85,	627.85,	677.85,	727.85,	777.85,	827.85,	877.85,	927.85,	977.85,	1027.85, 1077.85, 1084.85 ];
KC = [388.23 ,386.52 ,385.47, 384.59, 383.6, 382.34, 380.69, 378.59, 376.06, 373.12, 369.85, 366.33, 362.67, 358.98, 355.33, 351.8,	348.44,	345.26,	342.2, 339.16, 335.98, 332.39, 328.07,327.39];
CCK = polyfit( TCK, KC, 4);
kcf = int ( CCK(1)*T^4 + CCK(2)*T^3 + CCK(3)*T^2 + CCK(4)*T + CCK(5));

% Resistivity Copper 102 
TCP = [0, 19.85, 26.85,	76.85, 126.85, 226.85, 326.85, 426.85, 526.85, 626.85, 726.85, 826.85, 926.85, 1026.85, 1084.45 ];
PC = [1.5430E-08	1.6780E-08	1.7250E-08	2.0630E-08	2.4020E-08	3.0900E-08	3.7920E-08	4.5140E-08	5.2620E-08	6.0410E-08	6.8580E-08	7.7170E-08	8.6260E-08	9.5920E-08	1.0171E-07];
CCP = polyfit( TCP, PC, 4);
pcf = int ( CCP(1)*T^4 + CCP(2)*T^3 + CCP(3)*T^2 + CCP(4)*T + CCP(5));

% Conductivity Al2O3
TAK = [19.85 ,37.8298, 55.8096,	73.78939, 91.76919,	109.749, 127.7288, 145.7086, 163.6884, 181.6682, 199.648, 217.6278, 235.6076, 253.5874,	271.5672, 289.547, 307.5268, 325.5066, 343.4864, 361.4662 ...
    379.446, 397.4258, 415.4056, 433.3854, 451.3652, 469.3449, 487.3247, 505.3045, 523.2843, 541.2641, 559.2439,577.2237, 595.2035, 613.1833, 631.1631, 649.1429, 667.1227, 685.1025, 703.0823,	721.0621 ...
    739.0419, 757.0217, 775.0015, 792.9813,	810.9611, 828.9409,	846.9207, 864.9005,	882.8803, 900.8601,	918.8399, 936.8197,	954.7995,	972.7793,	990.7591,	1008.739,	1026.719,	1044.698,	1062.678 ...
    1080.658	1098.638	1116.618	1134.597	1152.577	1170.557	1188.537	1206.517	1224.496	1242.476	1260.456	1278.436	1296.416	1314.395	1332.375	1350.355	1368.335 ...
    1386.315	1404.294	1422.274	1440.254	1458.234	1476.214	1494.193	1512.173	1530.153	1548.133	1566.113	1584.092	1602.072	1620.052	1638.032	1656.012	1673.991 ...
    1691.971	1709.951	1727.931	1745.911	1763.89	1781.87	1799.85];

KA = [35.4396	33.51727	31.9325	30.42353	28.98777	27.6227	26.32582	25.0947	23.92694	22.82018	21.77212	20.78048	19.84306	18.95768	18.1222	17.33454	16.59267	15.89457	15.23832 ...
    14.62198	14.04371	13.50169	12.99413	12.51932	12.07557	11.66125	11.27474	10.91452	10.57906	10.26692	9.976668	9.70694	9.456409	9.223791	9.007849	8.80739	8.621267	...
    8.448375	8.287655	8.138092	7.998718	7.868607	7.746878	7.632696	7.525269	7.42385	7.327738	7.236275	7.148849	7.064891	6.983877	6.90533	6.828814	6.75394	6.680364	6.607784 ...
    6.535945	6.464636	6.393691	6.322987	6.252448	6.18204	6.111777	6.041714	5.971953	5.902639	5.833964	5.766163	5.699515	5.634345	5.571022	5.509959	5.451616	5.396494 ...
    5.345142	5.298152	5.25616	5.219848	5.189943	5.167214	5.152478	5.146594	5.150467	5.165047	5.191326	5.230344	5.283184	5.350973	5.434884	5.536134	5.655984	5.795741 ...
    5.956756	6.140425	6.348187	6.581527	6.841975	7.131104	7.450535	7.801928 ];
CAK = polyfit( TAK, KA, 4);
kaf = int ( CAK(1)*T^4 + CAK(2)*T^3 + CAK(3)*T^2 + CAK(4)*T + CAK(5));
 
 
syms kn kp sn sp pp pn T % k: Conductivity, s: Seebeck, p: resistivity, p: p-leg, n:n-leg

% Material 1
kn = symfun( knf, T); % Temperature dependent k for n-leg
kp = symfun (kpf, T); % Temperature dependent k for p-leg
sn = symfun (snf, T); % Temperature dependent s for n-leg
sp = symfun (spf, T); % Temperature dependent k for p-leg
pp = symfun (ppf, T); % Temperature dependent k for p-leg
pn = symfun (pnf, T); % Temperature dependent k for p-leg
% New code
kc = symfun (kcf, T); % Temperature dependent k for Copper 102
pc = symfun (pcf, T); % Temperature dependent p for Copper 102
ka = symfun (kaf, T); % Temperature dependent k for Al2O3
% New code
thn = symfun( thnf, T); % Temeperature dependent t for n-leg
thp = symfun( thpf, T); % Temeperature dependent t for n-leg

% Material 2
kn2 = symfun( knf2, T); % Temperature dependent k for n-leg
kp2 = symfun (kpf2, T); % Temperature dependent k for p-leg
sn2 = symfun (snf2, T); % Temperature dependent s for n-leg
sp2 = symfun (spf2, T); % Temperature dependent s for p-leg
pp2 = symfun (ppf2, T); % Temperature dependent p for p-leg
pn2 = symfun (pnf2, T); % Temperature dependent p for p-leg
% New code
kc = symfun (kcf, T); % Temperature dependent k for Copper 102
pc = symfun (pcf, T); % Temperature dependent p for Copper 102
ka = symfun (kaf, T); % Temperature dependent k for Al2O3

% Integral Averages 

knm = (kn(Th) - kn (Tc))/(Th-Tc);
knm = double(knm);
kpm = (kp(Th) - kp (Tc))/(Th-Tc);
kpm = double(kpm);
snm =  (sn(Th) - sn(Tc))/(Th-Tc);
snm = double(snm);
spm =  (sp(Th) - sp(Tc))/(Th-Tc);
spm = double(spm);
ppm =  (pp(Th) - pp(Tc))/(Th-Tc);
ppm = double(ppm);
pnm =  (pn(Th) - pn(Tc))/(Th-Tc);
pnm = double(pnm);

Th2 = TNS2(end);
Tc2 = TNS2(1);
knm2 = (kn2(Th2) - kn2 (Tc2))/(Th2-Tc2);
knm2 = double(knm2);
kpm2 = (kp2(Th2) - kp2 (Tc2))/(Th2-Tc2);
kpm2 = double(kpm2);
snm2 =  (sn2(Th2) - sn2(Tc2))/(Th2-Tc2);
snm2 = double(snm2);
spm2 =  (sp2(Th2) - sp2(Tc2))/(Th2-Tc2);
spm2 = double(spm2);
ppm2 =  (pp2(Th2) - pp2(Tc2))/(Th2-Tc2);
ppm2 = double(ppm2);
pnm2 =  (pn2(Th2) - pn2(Tc2))/(Th2-Tc2);
pnm2 = double(pnm2);


% New code
kcm =  (kc(Th) - kc(Tc))/(Th-Tc);
kcm = double(kcm);
pcm =  (pc(Th) - pc(Tc))/(Th-Tc);
pcm = double(pcm);
%
kam =  (ka(Th) - ka(Tc))/(Th-Tc);
kam = double(kam);



rc1 = (t3t/(kcm*t3l*t3w/2));
ra1 = (t2t/(kam*t2l*t2w/2));
rc2 = (t1t/(kcm*t1l*t1w/2));
cp =  (kpm* lp*tp)/hp;  
rthp = 1/cp;
rthhp2 = (hp2)/(kpm2*lp2*tp2);
cn =  (knm* ln*tn)/hn; 
rthn = 1/cn;
rthhn2 = (hn2)/(knm2*ln2*tn2);

ra2 = (b2t/(kam.*b2l*b2w/2));
rc3 = (b1t/(kcm*b1l*b1w/2));
rc4 = (b3t/(kcm*b3l*b3w/2));

Rtot = rc1 + ra1 + rc2 + rthp + rthhp2+ rc3 + ra2 + rc4;

P1 = rc1/Rtot;
P2 = ra1/Rtot;
P3 = rc2/Rtot;
P4 = rthp/Rtot;
P5 = rthhp2/Rtot;
P6 = rc3/Rtot;
P7 = ra2/Rtot;
P8 = rc4/Rtot;

DelTP1 =  (P1*(Th-Tc));
DelTP2 = (P2*(Th-Tc));
DelTP3 = (P3*(Th-Tc));
DelTP4 = (P4*(Th-Tc));
DelTP5 = (P5*(Th-Tc));
DelTP6 = (P6*(Th-Tc));
DelTP7 = (P7*(Th-Tc));
DelTP8 = (P8*(Th-Tc));

DeltT = DelTP1+ DelTP2 + DelTP3+ DelTP4+ DelTP5+ DelTP6 + DelTP7 + DelTP8;

Tcu1 = Th - DelTP1;
Ta1 = Tcu1 - DelTP2;
Tcu2 = Ta1 - DelTP3;
Tpleg = Tcu2 - DelTP4;
Tpleg2 = Tpleg - DelTP5;
Tcu3 = Tpleg2 - DelTP6;
Ta2 = Tcu3 - DelTP7;
Tcu4 = Ta2 - DelTP8;

Tp = [ Th Tcu1 Ta1 Tcu2 Tpleg Tpleg2 Tcu3 Ta2  Tc];
%
Rtotn = rc1 + ra1 + rc2 + rthn + rthhn2 + rc3 + ra2 + rc4;

N1 = rc1/Rtotn;
N2 = ra1/Rtotn;
N3 = rc2/Rtotn;
N4 = rthn/Rtotn;
N5 = rthhn2/Rtotn;
N6 = rc3/Rtotn;
N7 = ra2/Rtotn;
N8 = rc4/Rtotn;

DelTN1 =  (N1*(Th-Tc));
DelTN2 = (N2*(Th-Tc));
DelTN3 = (N3*(Th-Tc));
DelTN4 = (N4*(Th-Tc));
DelTN5 = (N5*(Th-Tc));
DelTN6 = (N6*(Th-Tc));
DelTN7 = (N7*(Th-Tc));
DelTN8 = (N8*(Th-Tc));


Tcu1n = Th - DelTN1;
Ta1n = Tcu1n - DelTN2;
Tcu2n = Ta1n - DelTN3;
Tnleg = Tcu2n - DelTN4;
Tnleg2 = Tnleg - DelTN5;
Tcu3n = Tnleg2 - DelTN6;
Ta2n = Tcu3n - DelTN7;
Tcu4n = Ta2n - DelTN8;

Tn = [ Th Tcu1n Ta1n Tcu2n Tnleg Tcu3n Ta2n  Tc];

TotalHeight = t3t + hp +hp2 + t2t + t1t + b3t + b2t + b1t;
%
t3telems = round((t3t/TotalHeight)*elems);
t2telems = round((t2t/TotalHeight)*elems);
t1telems = round((t1t/TotalHeight)*elems);
hpelems = round((hp/TotalHeight)*elems);
hp2elems = round((hp2/TotalHeight)*elems);
b3telems = round((b3t/TotalHeight)*elems);
b2telems = round((b2t/TotalHeight)*elems);
b1telems = round((b1t/TotalHeight)*elems);
summedelems = sum([t3telems , t2telems, t1telems, hpelems, hp2elems, b3telems, b2telems ,b1telems]);

while ( summedelems ~= elems)
%     fprintf ('Problem with elems\n');
    if (summedelems > elems & t3telems > b3telems)
        t3telems = t3telems -1;
        summedelems = sum([t3telems , t2telems, t1telems, hpelems, b3telems, b2telems ,b1telems]);
    elseif (summedelems > elems)
        b3telems = b3telems -1;
        summedelems = sum([t3telems , t2telems, t1telems, hpelems, b3telems, b2telems ,b1telems]);
    elseif (summedelems < elems & t3telems > b3telems)
        b3telems = b3telems +1;
        summedelems = sum([t3telems , t2telems, t1telems, hpelems, b3telems, b2telems ,b1telems]);
    else
        t3telems = t3telems +1;
        summedelems = sum([t3telems , t2telems, t1telems, hpelems, b3telems, b2telems ,b1telems]);

    end
end


dt = (Th - Tcu1)/t3telems;
dtn = (Th - Tcu1n)/t3telems;
%
for i = 1: t3telems 
    
    T_t3t(i+1) = Th - dt*(i);
    T_t3tn(i+1) = Th - dtn*(i);
    
end
T_t3t(1) = Th;
T_t3tn(1) = Th;

dt = (Tcu1 - Ta1)/t2telems;
dtn = (Tcu1n - Ta1n)/t2telems;


for i = 1: t2telems 
    
    T_t2t(i+1) = Tcu1 - dt*(i);
    T_t2tn(i+1) = Tcu1n - dtn*(i);
    
end
T_t2t(1) = Tcu1;
T_t2tn(1) = Tcu1n;

dt = (Ta1 - Tcu2)/t1telems;
dtn = (Ta1n - Tcu2n)/t1telems;


for i = 1: t1telems 
    
    T_t1t(i+1) = Tcu2 - dt*(i);
    T_t1tn(i+1) = Tcu2n - dt*(i);

    
end
T_t1t(1) = Ta1;
T_t1tn(1) = Ta1n;

dt = (Tcu2 - Tpleg)/hpelems;
dtn = (Tcu2n - Tnleg)/hpelems;


for i = 1: hpelems 
    
    T_hp(i+1) = Tcu2 - dt*(i);
    T_hn(i+1) = Tcu2n - dtn*(i);
    
end

T_hp(1) = Tcu2;
T_hn(1) = Tcu2n;

dt = (Tpleg - Tpleg2)/hp2elems;
dtn = (Tnleg - Tnleg2)/hp2elems;
% Tp = [ Th Tcu1 Ta1 Tcu2 Tpleg Tpleg2 Tcu3 Ta2  Tc]


for i = 1: hp2elems 
    
    T_hp2(i+1) = Tpleg - dt*(i);
    T_hn2(i+1) = Tnleg - dtn*(i);
    
end

T_hp2(1) = Tpleg;
T_hn2(1) = Tnleg;

dt = (Tpleg2 - Tcu3)/b1telems;
dtn = (Tnleg2 - Tcu3n)/b1telems;

for i = 1: b1telems 
     
    T_b1t(i+1) = Tpleg2 - dt*(i);
    T_b1tn(i+1) = Tnleg2 - dtn*(i);
    
end
T_b1t(1) = Tpleg2;
T_b1tn(1) = Tnleg2;

dt = (Tcu3 - Ta2)/b2telems;
dtn = (Tcu3n - Ta2n)/b2telems;

for i = 1: b2telems 
    
    T_b2t(i+1) = Tcu3 - dt*(i);
    T_b2tn(i+1) = Tcu3n - dt*(i);

    
end
T_b2t(1) = Tcu3;
T_b2tn(1) = Tcu3n;

dt = (Ta2 - Tc)/b3telems;
dtn = (Ta2n - Tc)/b3telems;



for i = 1: b3telems 
    
    T_b3t(i+1) = Ta2 - dt*(i);
    T_b3tn(i+1) = Ta2 - dt*(i);

    
end
T_b3t(1) = Ta2;
T_b3tn(1) = Ta2n;


Tp = [ T_t3t T_t2t T_t1t T_hp T_hp2 T_b1t T_b2t T_b3t ];
Tn = [ T_t3tn T_t2tn T_t1tn T_hn T_hn2 T_b1tn T_b2tn T_b3tn ];

 
% Energy Generation

ppm =  (pp(T_hp(1)) - pp(T_hp(end)))/(T_hp(1)-T_hp(end));
ppm = double(ppm);
pnm =  (pn(T_hn(1)) - pn(T_hn(end)))/(T_hn(1)-T_hn(end));
pnm = double(pnm);

ppm2 =  (pp2(T_hp2(1)) - pp2(T_hp2(end)))/(T_hp2(1)-T_hp2(end));
ppm2 = double(ppm2);
pnm2 =  (pn2(T_hn2(1)) - pn2(T_hn2(end)))/(T_hn2(1)-T_hn2(end));
pnm2 = double(pnm2);

rp = (ppm * hp )/(lp*tp); % [Ohms] Resistance p-leg
rn = (pnm * hn )/(ln*tn); % [Ohms] Resistance n-leg

rp2 = (ppm2 * hp2 )/(lp2*tp2); % [Ohms] Resistance p-leg
rn2 = (pnm2 * hn2 )/(ln2*tn2); % [Ohms] Resistance n-leg

pt1tm =  (pc(T_t1t(1)) - pc(T_t1t(end)))/(T_t1t(1) - T_t1t(end) );
pt1tm = double(pt1tm);

pb1tpm =  (pc(T_b1t(1)) - pc(T_b1t(end)))/(T_b1t(1) - T_b1t(end) );
pb1tpm = double(pb1tpm);

pb1tnm =  (pc(T_b1tn(1)) - pc(T_b1tn(end)))/(T_b1tn(1) - T_b1tn(end) );
pb1tnm = double(pb1tnm);



snm =  (sn(T_hn(1)) - sn(T_hn(end)))/(T_hn(1)-T_hn(end));
snm = double(snm);
spm =  (sp(T_hp(1)) - sp(T_hp(end)))/(T_hp(1)-T_hp(end));
spm = double(spm);

snm2 =  (sn2(T_hn2(1)) - sn2(T_hn2(end)))/(T_hn2(1)-T_hn2(end));
snm2 = double(snm2);
spm2 =  (sp2(T_hp2(1)) - sp2(T_hp2(end)))/(T_hp2(1)-T_hp2(end));
spm2 = double(spm2);
% Contact Resistance

CR_n = C_pho/(tp*100)*(lp*100) + 2*C_pho/(tp2*100)*(lp2*100) ;
CR_p = C_pho/(tp*100)*(lp*100) + 2*C_pho/(tn2*100)*(ln2*100) ;
CR = (CR_n + CR_p);

% 

Rt1t = (pt1tm* t1w)/(t1l*t1t);
Rb1tp = (pb1tpm* b1w)/(b1l*b1t);
Rb1tn =  (pb1tnm* b1w)/(b1l*b1t) ;
Rt = rp +rn + rp2 + rn2  + Rt1t + Rb1tp +Rb1tn + CR ;% [Ohms] Resistance Total of legs
RL = Rt ; % [Ohms] Load resistance is set to internal resistance
R = RL + Rt ; % [Ohms] 
S =  spm +spm2  - snm -snm2 ; % Seebeck of couple
I = (S*(T_hp(1) -T_hp2(end)))/(R);% [A] Current through the unicouple

%  I = 3.64
% 

t1te_h = t1t/t1telems;

Jh_t1t = (I^2)* Rt1t;
Jh_b1tp = (I^2)*Rb1tp;
Jh_b1tn = (I^2)*Rb1tn;

Jh_t1t_e = Jh_t1t/(t1telems*2);
Jh_b1tp_e = Jh_t1t/(b1telems);
Jh_b1tn_e = Jh_t1t/(b1telems);

shp = hp/hpelems;
shn = hn/hpelems;


for j = 1:hpelems
  
    knm(j) = (kn(T_hn(j)) - kn(T_hn(j+1)))/(T_hn(j)-T_hn(j+1));
    kpm(j) = (kp(T_hp(j)) - kp(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));

    snm(j) =  (sn(T_hn(j)) - sn(T_hn(j+1)))/(T_hn(j)-T_hn(j+1));
    spm(j) =  (sp(T_hp(j)) - sp(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));
    ppm(j) =  (pp(T_hp(j)) - pp(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));
    pnm(j) =  (pn(T_hn(j)) - pn(T_hn(j+1)))/(T_hn(j)-T_hn(j+1));
    thpm(j) =  (thp(T_hp(j)) - thp(T_hp(j+1)))/(T_hp(j)-T_hp(j+1)); % new code
    thpm(j) = double(thpm(j));
    thnm(j) =  (thn(T_hn(j)) - thn(T_hn(j+1)))/(T_hn(j)-T_hn(j+1)); % new code
    thnm(j) = double(thnm(j));

    rp1(j) = (ppm(j)* shp )/(lp*tp);  % [Ohms] Resistance p-leg
    rn1(j) = (pnm(j)* shn )/(ln*tn); % [Ohms] Resistance n-leg
    cp(j) =  (kpm(j)* lp*tp)/shp; %  Thermal Conductance of p-leg
    cn(j) = (knm(j)* ln*tn)/shn; % Thermal Conductance of n-leg
    
    Qhp(j) = spm(j)*T_hp(j)*I + cp(j)*(T_hp(j)-T_hp(j+1)) -(0.5)*(I^2)*rp1(j); % - 0.5*thpm(j)*I;  % Heat input to hot side
    Qcp(j) = spm(j)*T_hp(j+1)*I + cp(j)*(T_hp(j)-T_hp(j+1)) + 0.5*(I^2)*rp1(j); % + 0.5*thpm(j)*I; % Heat rejection from cold side
    Pp(j) = double((Qhp(j)- Qcp(j))) ;% Power output per leg
   
    Qhn(j) = abs(snm(j))*T_hn(j)*I + cn(j)*(T_hn(j)-T_hn(j+1)) - (0.5)*(I^2)*rn1(j) ; %- 0.5*(thnm(j))*I ; % Heat input to hot side
    Qcn(j) = abs(snm(j))*T_hn(j+1)*I + cn(j)*(T_hn(j)-T_hn(j+1))  + 0.5*(I^2)*rn1(j) ; %+ 0.5*((thnm(j)))*I ;% Heat rejection from cold side
    Pn(j) = double(Qhn(j)- Qcn(j)); % Power output per leg


end

InitIt_Powp = sum(Pp); 
InitIt_Pown = sum(Pn); 

shp = hp2/hp2elems;
shn = hn2/hp2elems;

for j = 1:hp2elems
  
    knm(j) = (kn2(T_hn2(j)) - kn2(T_hn2(j+1)))/(T_hn2(j)-T_hn2(j+1));
    kpm(j) = (kp2(T_hp2(j)) - kp2(T_hp2(j+1)))/(T_hp2(j)-T_hp2(j+1));

    snm(j) =  (sn2(T_hn2(j)) - sn2(T_hn2(j+1)))/(T_hn2(j)-T_hn2(j+1));
    spm(j) =  (sp2(T_hp2(j)) - sp2(T_hp2(j+1)))/(T_hp2(j)-T_hp2(j+1));
    ppm(j) =  (pp2(T_hp2(j)) - pp2(T_hp2(j+1)))/(T_hp2(j)-T_hp2(j+1));
    pnm(j) =  (pn2(T_hn2(j)) - pn2(T_hn2(j+1)))/(T_hn2(j)-T_hn2(j+1));
%     thpm(j) =  (thp(T_hp(j)) - thp(T_hp(j+1)))/(T_hp(j)-T_hp(j+1)); % new code
%     thpm(j) = double(thpm(j));
%     thnm(j) =  (thn(T_hn(j)) - thn(T_hn(j+1)))/(T_hn(j)-T_hn(j+1)); % new code
%     thnm(j) = double(thnm(j));

    rp1(j) = (ppm(j)* shp )/(lp*tp);  % [Ohms] Resistance p-leg
    rn1(j) = (pnm(j)* shn )/(ln*tn); % [Ohms] Resistance n-leg
    cp(j) =  (kpm(j)* lp*tp)/shp; %  Thermal Conductance of p-leg
    cn(j) = (knm(j)* ln*tn)/shn; % Thermal Conductance of n-leg
    
    Qhp(j) = spm(j)*T_hp(j)*I + cp(j)*(T_hp(j)-T_hp(j+1)) -(0.5)*(I^2)*rp1(j); % - 0.5*thpm(j)*I;  % Heat input to hot side
    Qcp(j) = spm(j)*T_hp(j+1)*I + cp(j)*(T_hp(j)-T_hp(j+1))+ 0.5*(I^2)*rp1(j); % + 0.5*thpm(j)*I; % Heat rejection from cold side
    Pp2(j) = double((Qhp(j)- Qcp(j))) ;% Power output per leg
   
    Qhn(j) = abs(snm(j))*T_hn(j)*I + cn(j)*(T_hn(j)-T_hn(j+1)) - (0.5)*(I^2)*rn1(j) ; %- 0.5*(thnm(j))*I ; % Heat input to hot side
    Qcn(j) = abs(snm(j))*T_hn(j+1)*I + cn(j)*(T_hn(j)-T_hn(j+1)) + 0.5*(I^2)*rn1(j) ; %+ 0.5*((thnm(j)))*I ;% Heat rejection from cold side
    Pn2(j) = double(Qhn(j)- Qcn(j)); % Power output per leg


end

InitIt_Powp2 = sum(Pp2); 
InitIt_Pown2 = sum(Pn2);



for i = 1:t3telems %  t3t (Very Top Cu Layer) Copper Conductivity and Area
    
    k_t3tp_e(i) = double((kc(T_t3t(i)) - kc(T_t3t(i+1)))/(T_t3t(i) - T_t3t(i+1)));
    k_t3tn_e(i) = double((kc(T_t3tn(i)) - kc(T_t3tn(i+1)))/(T_t3tn(i) - T_t3tn(i+1)));
    A_t3t(i) = t3w*t3l;

end

for i = 1:t2telems % t22 Top Alumina Conductivity and Area
    
    k_t2tp_e(i) = double((ka(T_t2t(i)) - ka(T_t2t(i+1)))/(T_t2t(i) - T_t2t(i+1)));
    k_t2tn_e(i) = double((ka(T_t2tn(i)) - ka(T_t2tn(i+1)))/(T_t2tn(i) - T_t2tn(i+1)));
    A_t2tp(i) = 0.5*(t2w*t2l);
    A_t2tn(i) = 0.5*(t2w*t2l);

end

for i = 1:t1telems % t1t ( Copper Layer Connecting Legs ) Conductivity and Area
    
    k_t1tp_e(i) = double((kc(T_t1t(i)) - kc(T_t1t(i+1)))/(T_t1t(i) - T_t1t(i+1)));
    k_t1tn_e(i) = double((kc(T_t1tn(i)) - kc(T_t1tn(i+1)))/(T_t1tn(i) - T_t1tn(i+1)));
    A_t1tp(i) = 0.5*(t1w*t1l);
    A_t1tn(i) = 0.5*(t1w*t1l);
end

for i = 1:hpelems % Leg conductivity and area
    
    k_hp_e(i) = double((kp(T_hp(i)) - kp(T_hp(i+1)))/(T_hp(i) - T_hp(i+1)));
    k_hn_e(i) = double((kn(T_hn(i)) - kn(T_hn(i+1)))/(T_hn(i) - T_hn(i+1)));
    A_p(i) = lp*tp;
    A_n(i) = ln*tn;
    
    
end

for i = 1:hp2elems % Leg conductivity and area
    
    k_hp2_e(i) = double((kp2(T_hp2(i)) - kp2(T_hp2(i+1)))/(T_hp2(i) - T_hp2(i+1)));
    k_hn2_e(i) = double((kn2(T_hn2(i)) - kn2(T_hn2(i+1)))/(T_hn2(i) - T_hn2(i+1)));
    A_p2(i) = lp2*tp2;
    A_n2(i) = ln2*tn2;
    
    
end

for i = 1:b1telems % Bottom Copper Conductivity and Area
    
    k_b1tp_e(i) = double((kc(T_b1t(i)) - kc(T_b1t(i+1)))/(T_b1t(i) - T_b1t(i+1)));
    k_b1tn_e(i) = double((kc(T_b1tn(i)) - kc(T_b1tn(i+1)))/(T_b1tn(i) - T_b1tn(i+1)));
    A_b1t(i) = b1w*b1l;
end

for i = 1:b2telems % Bottom Alumina Conductivity and Area
    
    k_b2tp_e(i) = double((ka(T_b2t(i)) - ka(T_b2t(i+1)))/(T_b2t(i) - T_b2t(i+1)));
    k_b2tn_e(i) = double((ka(T_b2tn(i)) - ka(T_b2tn(i+1)))/(T_b2tn(i) - T_b2tn(i+1)));
    A_b2tp(i) = 0.5*(b2w*b2l);
    A_b2tn(i) = 0.5*(b2w*b2l);
end


for i = 1:b3telems % Bottomost Copper Conductivity and Area
    
    k_b3tp_e(i) = double((kc(T_b3t(i)) - kc(T_b3t(i+1)))/(T_b3t(i) - T_b3t(i+1)));
    k_b3tn_e(i) = double((kc(T_b3tn(i)) - kc(T_b3tn(i+1)))/(T_b3tn(i) - T_b3tn(i+1)));
    A_b3tp(i) = 0.5*(b3w*b3l);
    A_b3tn(i) = 0.5*(b3w*b3l);

end

kp_elems = [ k_t3tp_e, k_t2tp_e, k_t1tp_e, k_hp_e, k_hp2_e, k_b1tp_e, k_b2tp_e, k_b3tp_e];
kn_elems = [ k_t3tn_e, k_t2tn_e, k_t1tn_e, k_hn_e, k_hn2_e, k_b1tn_e, k_b2tn_e, k_b3tn_e];

Ap_elems = [ A_t3t, A_t2tp, A_t1tp, A_p, A_p2,  A_b1t, A_b2tp, A_b3tp];
An_elems = [ A_t3t, A_t2tn, A_t1tn, A_n, A_n2, A_b1t, A_b2tn, A_b3tn];

% ke = (A*k)/(6*l)*[14, -16, 2; -16, 32, -16; 2, -16, 14];
Kp = sparse(numn, numn);
Kn = sparse(numn, numn);

el_l = TotalHeight/elems;

% K - Matrix
for i=1:2:numn-2
    
    if (i == 1)
        j = i ;
    else
        j =i-(prvj);
    end
    prvj = j;
    
    ke_p = (kp_elems(j))*(Ap_elems(j))/(6*el_l)*[14, -16, 2; -16, 32, -16; 2, -16, 14];
    ke_n = (kn_elems(j))*(An_elems(j))/(6*el_l)*[14, -16, 2; -16, 32, -16; 2, -16, 14];
    
    dof = [ i, i+1, i+2];
    Kp(dof, dof) = ke_p +  Kp(dof, dof);
    Kn(dof, dof) = ke_n +  Kn(dof, dof);

end
% 


% F - Vector
C0 = t3telems;
C1 = t3telems + t2telems; % Ending Element for t2t
C2 = C1 + t1telems; % Ending Element for t1t
C3 = C2 + hpelems; % Ending element for legs
C4 = C3 + hp2elems;% Ending element for material 2
C5 = C4 + b1telems; % Ending element for b1t
C6 = C5 + b2telems;

for i = C1+1: C2
    
    E_gen_t1t(i) = Jh_t1t_e/(Ap_elems(i) * el_l);
    
end

for i = 1:hpelems

    E_genp(i) = -Pp(i)./(Ap_elems(i+C2).*el_l);
    E_genn(i) = -Pn(i)/(Ap_elems(i+C2)*el_l);

end

for i = 1:hp2elems

    E_genp2(i) = -Pp2(i)./(Ap_elems(i+C2).*el_l);
    E_genn2(i) = -Pn2(i)/(Ap_elems(i+C2)*el_l);

end
%  E_gen_b1t = zeros(b1telems +b2telems + b3telems,1)';
for i = C4+1: C5
    
    j= i -C4;
    E_gen_b1tp(j) = Jh_b1tp_e./(Ap_elems(i)* el_l);
    E_gen_b1tn(j) = Jh_b1tn_e./(Ap_elems(i)* el_l);

end



E_gen_b2t = zeros(b2telems,1)';
E_gen_b3t = zeros(b3telems,1)';

E_Gen_pL = [ E_gen_t1t, E_genp, E_genp2, E_gen_b1tp, E_gen_b2t, E_gen_b3t];
E_Gen_nL = [ E_gen_t1t, E_genn, E_genn2, E_gen_b1tn, E_gen_b2t, E_gen_b3t];
    
prvj =1;
for i = 3:2: numn-2
    
    j = i - prvj;
    prvj =j;
    fe_p = (E_Gen_pL(j)*Ap_elems(j)*el_l/6)*[1, 4, 1]';
    fe_n = (E_Gen_nL(j)*An_elems(j)*el_l/6)*[1, 4, 1]';

    dof = [i, i+1, i+2];
    F_p(dof) =fe_p;
    F_p(i) = fe_p(1) + fe_p(3);
           
    F_n(dof) =fe_n;
    F_n(i) = fe_n(1) + fe_n(3);
end

F_p(1) = Th;
F_n(1) = Th;
F_p(2) = F_p(2) - Kp(2,1)*Th;
F_n(2) = F_n(2) - Kn(2,1)*Th;
F_p(3) = F_p(3) - Kp(3,1)*Th;
F_n(3) = F_n(3) - Kn(3,1)*Th;

F_p(end-2) = F_p(end -2) - Kp(end-2,end)*Tc;
F_p(end -1) = F_p(end -1) - Kp(end-1,end)*Tc;
F_p(end) = Tc;

F_n(end-2) = F_n(end -2) - Kn(end-2,end)*Tc;
F_n(end -1) = F_n(end -1) - Kn(end-1,end)*Tc;
F_n(end) = Tc;

for i = 1: numn
    Kp(numn,i) = 0;
    Kp(i,numn) = 0;
    Kp(numn,numn) =1;
    Kp(1,i) = 0;
    Kp(i,1) = 0;
    
    Kn(numn,i) = 0;
    Kn(i,numn) = 0;
    Kn(numn,numn) =1;
    Kn(1,i) = 0;
    Kn(i,1) = 0;

end
Kp(1,1) =1;
Kn(1,1) =1;


Tempr_p = Kp\F_p';
Tempr_n = Kn\F_n';

%

% C0 = t3telems;
% C1 = t3telems + t2telems; % Ending Element for t2t
% C2 = C1 + t1telems; % Ending Element for t1t
% C3 = C2 + hpelems; % Ending element for legs
% C4 = C3 + hp2elems;% Ending element for material 2
% C5 = C4 + b1telems; % Ending element for b1t
% C6 = C5 + b2telems;


ppm =  (pp(Tempr_p(C2*2 +1)) - pp(Tempr_p(C3*2 +1)))/(Tempr_p(C2*2 +1)-Tempr_p(C3*2 +1));
ppm = double(ppm);
pnm =  (pn(Tempr_p(C2*2 +1)) - pn(Tempr_p(C3*2 +1)))/(Tempr_p(C2*2 +1)-Tempr_p(C3*2 +1));
pnm = double(pnm);


ppm2 =  (pp2(Tempr_p(C3*2 +1)) - pp2(Tempr_p(C4*2 +1)))/(Tempr_p(C3*2 +1)-Tempr_p(C4*2 +1));
ppm2 = double(ppm2);
pnm2 =  (pn2(Tempr_p(C3*2 +1)) - pn2(Tempr_p(C4*2 +1)))/(Tempr_p(C3*2 +1)-Tempr_p(C4*2 +1));
pnm2 = double(pnm2);

pt1tm =  (pc(Tempr_p(C1*2 +1 )) - pc(Tempr_p(C2*2 +1)))/(Tempr_p(C1*2 +1) - Tempr_p(C2*2 +1) );
pt1tm = double(pt1tm);

pb1tpm =  (pc(Tempr_p(C3*2 +1)) - pc(Tempr_p(C4*2 +1)))/(Tempr_p(C3*2 +1) - Tempr_p(C4*2 +1) );
pb1tpm = double(pb1tpm);

pb1tnm =  (pc(Tempr_n(C3*2 +1)) - pc(Tempr_n(C4*2 +1)))/(Tempr_n(C3*2 +1) - Tempr_n(C4*2 +1) );
pb1tnm = double(pb1tnm);

snm =  (sn(Tempr_n(C2*2 +1)) - sn(Tempr_n(C3*2 +1)))/(Tempr_n(C2*2 +1)-Tempr_n(C3*2 +1));
snm = double(snm);
spm =  (sp(Tempr_p(C2*2 +1)) - sp(Tempr_p(C3*2 +1)))/(Tempr_p(C2*2 +1)-Tempr_p(C3*2 +1));
spm = double(spm);

snm2 =  (sn2(Tempr_n(C3*2 +1)) - sn2(Tempr_n(C4*2 +1)))/(Tempr_n(C3*2 +1)-Tempr_n(C4*2 +1));
snm2 = double(snm2);
spm2 = (sp2(Tempr_p(C3*2 +1)) - sp2(Tempr_p(C4*2 +1)))/(Tempr_p(C3*2 +1)-Tempr_p(C4*2 +1));
spm2 = double(spm2);

rp = (ppm * hp )/(lp*tp); % [Ohms] Resistance p-leg
rn = (pnm * hn )/(ln*tn); % [Ohms] Resistance n-leg

rp2 = (ppm2 * hp2 )/(lp2*tp2); % [Ohms] Resistance p-leg
rn2 = (pnm2 * hn2 )/(ln2*tn2); % [Ohms] Resistance n-leg


Rt1t = (pt1tm* t1w)/(t1l*t1t);
Rb1tp = (pb1tpm* b1w)/(b1l*b1t);
Rb1tn =  (pb1tnm* b1w)/(b1l*b1t) ;
Rt = rp +rn + rp2+ rn2+ Rt1t + Rb1tp +Rb1tn + CR ;% [Ohms] Resistance Total of legs
RL = Rt ; % [Ohms] Load resistance is set to internal resistance
R = RL + Rt ; % [Ohms] 

S =  spm +spm2  -snm -snm2;  % Seebeck of couple
Voc_av = S*(Tempr_p(C2*2 +1) -Tempr_p(C4*2 +1));
I = S*(Tempr_p(C2*2 +1) -Tempr_p(C3*2 +1))./(R);% [A] Current through the unicouple

shp = hp/hpelems;
shn = hn/hpelems;

% I = 3.64;
%  I = 9.2876;
for k = 1:hpelems
  
    if ( k ==1)
        j = C2*2 + 1;
    else
        j = C2*2 +1 + (k-1)*2;
    end
    
    
    
    knm(k) = (kn(Tempr_p(j)) - kn(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    kpm(k) = (kp(Tempr_p(j)) - kp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));

    snm(k) =  (sn(Tempr_n(j)) - sn(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2));
    spm(k) =  (sp(Tempr_p(j)) - sp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    ppm(k) =  (pp(Tempr_p(j)) - pp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    pnm(k) =  (pn(Tempr_n(j)) - pn(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2));
%     thpm(k) =  (thp(Tempr_p(j)) - thp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2)); % new code
%     thpm(k) = double(thpm(k));
%     thnm(k) =  (thn(Tempr_n(j)) - thn(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2)); % new code
%     thnm(k) = double(thnm(k));

    rp1(k) = (ppm(k)* shp )/(lp*tp);  % [Ohms] Resistance p-leg
    rn1(k) = (pnm(k)* shn )/(ln*tn); % [Ohms] Resistance n-leg
    cp(k) =  (kpm(k)* lp*tp)/shp; %  Thermal Conductance of p-leg
    cn(k) = (knm(k)* ln*tn)/shn; % Thermal Conductance of n-leg
    
    V_oc_e(k) = (spm(k)-snm(k))*(Tempr_p(j) -Tempr_p(j+2));
    
    Qhp(k) = spm(k)*Tempr_p(j)*I  + cp(k)*(Tempr_p(j)-Tempr_p(j+2))-(0.5)*(I^2)*rp1(k) ; % - abs(0.5*thpm(k)*I);  % Heat input to hot side
    Qcp(k) = spm(k)*Tempr_p(j+2)*I + cp(k)*(Tempr_p(j)-Tempr_p(j+2))+ 0.5*(I^2)*rp1(k) ; %+ abs(0.5*thpm(k)*I); % Heat rejection from cold side
    Pp(k) = double((Qhp(k)- Qcp(k))) ;% Power output per leg
   
    Qhn(k) = abs(snm(k))*Tempr_p(j)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) - (0.5)*(I^2)*rn1(k); % - abs(0.5*(thnm(k))*I) ; % Heat input to hot side
    Qcn(k) = abs(snm(k))*Tempr_p(j+2)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) + 0.5*(I^2)*rn1(k); % + abs(0.5*((thnm(k)))*I) ;% Heat rejection from cold side
    Pn(k) = double(Qhn(k)- Qcn(k)); % Power output per leg


end

shp2 = hp2/hp2elems;
shn2 = hn2/hp2elems;

for k = 1:hp2elems
  
    if ( k ==1)
        j = C3*2 + 1;
    else
        j = C3*2 +1 + (k-1)*2;
    end
    
    
    
    knm(k) = (kn2(Tempr_p(j)) - kn2(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    kpm(k) = (kp2(Tempr_p(j)) - kp2(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));

    snm(k) =  (sn2(Tempr_n(j)) - sn2(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2));
    spm(k) =  (sp2(Tempr_p(j)) - sp2(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    ppm(k) =  (pp2(Tempr_p(j)) - pp2(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    pnm(k) =  (pn2(Tempr_n(j)) - pn2(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2));
    thpm(k) =  (thp(Tempr_p(j)) - thp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2)); % new code
    thpm(k) = double(thpm(k));
    thnm(k) =  (thn(Tempr_n(j)) - thn(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2)); % new code
    thnm(k) = double(thnm(k));

    rp12(k) = (ppm(k)* shp2 )/(lp*tp);  % [Ohms] Resistance p-leg
    rn12(k) = (pnm(k)* shn2 )/(ln*tn); % [Ohms] Resistance n-leg
    cp(k) =  (kpm(k)* lp*tp)/shp2; %  Thermal Conductance of p-leg
    cn(k) = (knm(k)* ln*tn)/shn2; % Thermal Conductance of n-leg
    
    V_oc2_e(k) = (spm(k)-snm(k))*(Tempr_p(j) -Tempr_p(j+2));
    
    Qhp2(k) = spm(k)*Tempr_p(j)*I  + cp(k)*(Tempr_p(j)-Tempr_p(j+2))-(0.5)*(I^2)*rp12(k) ; % - abs(0.5*thpm(k)*I);  % Heat input to hot side
    Qcp2(k) = spm(k)*Tempr_p(j+2)*I + cp(k)*(Tempr_p(j)-Tempr_p(j+2))+ 0.5*(I^2)*rp12(k) ; %+ abs(0.5*thpm(k)*I); % Heat rejection from cold side
    Pp2(k) = double((Qhp2(k)- Qcp2(k))) ;% Power output per leg
   
    Qhn2(k) = abs(snm(k))*Tempr_p(j)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) - (0.5)*(I^2)*rn12(k); % - abs(0.5*(thnm(k))*I) ; % Heat input to hot side
    Qcn2(k) = abs(snm(k))*Tempr_p(j+2)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) + 0.5*(I^2)*rn12(k); % + abs(0.5*((thnm(k)))*I) ;% Heat rejection from cold side
    Pn2(k) = double(Qhn2(k)- Qcn2(k)); % Power output per leg


end

Qhin = max(Qhp) + max(Qhn);
Rpl = sum(rp1);
Rnl = sum(rn1);
Rpl2 = sum(rp12);
Rnl2 = sum(rn12);

R_FE = Rpl + Rnl + Rpl2 + Rnl2 +  Rt1t + Rb1tn + Rb1tp +CR;
% TotalP = sum(Pp);
% TotalN = sum(Pn);
% TotalP2 = sum(Pp2);
% TotalN2 = sum(Pn2);
% Power = TotalP + TotalN + TotalP2 + TotalN2 - (Rt1t + Rb1tn + Rb1tp + CR)*I^2
V_oc = sum(V_oc_e) + sum(V_oc2_e);
I_voc = V_oc/(R_FE*2);
I = I_voc;
% Power_El = (V_oc^2)/(4*R_FE)
% Eff = Power_El*100/Qhin
% RLoad = linspace (0, 7.5*R_FE, 100);
% 
% VTEG = R_FE.*(V_oc)./(RLoad + R_FE);
% ITEG = V_oc./(RLoad + R_FE);
% ITEG(end +1) = 0;
% VTEG = ITEG*R_FE;
% Vdev = V_oc - (R_FE*ITEG);

shp = hp/hpelems;
shn = hn/hpelems;

for k = 1:hpelems
  
    if ( k ==1)
        j = C2*2 + 1;
    else
        j = C2*2 +1 + (k-1)*2;
    end
    
    
    
    knm(k) = (kn(Tempr_p(j)) - kn(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    kpm(k) = (kp(Tempr_p(j)) - kp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));

    snm(k) =  (sn(Tempr_n(j)) - sn(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2));
    spm(k) =  (sp(Tempr_p(j)) - sp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    ppm(k) =  (pp(Tempr_p(j)) - pp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    pnm(k) =  (pn(Tempr_n(j)) - pn(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2));
%     thpm(k) =  (thp(Tempr_p(j)) - thp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2)); % new code
%     thpm(k) = double(thpm(k));
%     thnm(k) =  (thn(Tempr_n(j)) - thn(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2)); % new code
%     thnm(k) = double(thnm(k));

    rp1(k) = (ppm(k)* shp )/(lp*tp);  % [Ohms] Resistance p-leg
    rn1(k) = (pnm(k)* shn )/(ln*tn); % [Ohms] Resistance n-leg
    cp(k) =  (kpm(k)* lp*tp)/shp; %  Thermal Conductance of p-leg
    cn(k) = (knm(k)* ln*tn)/shn; % Thermal Conductance of n-leg
    
    V_oc_e(k) = (spm(k)-snm(k))*(Tempr_p(j) -Tempr_p(j+2));
    
    Qhp(k) = spm(k)*Tempr_p(j)*I  + cp(k)*(Tempr_p(j)-Tempr_p(j+2))-(0.5)*(I^2)*rp1(k) ; % - abs(0.5*thpm(k)*I);  % Heat input to hot side
    Qcp(k) = spm(k)*Tempr_p(j+2)*I + cp(k)*(Tempr_p(j)-Tempr_p(j+2))+ 0.5*(I^2)*rp1(k) ; %+ abs(0.5*thpm(k)*I); % Heat rejection from cold side
    Pp(k) = double((Qhp(k)- Qcp(k))) ;% Power output per leg
   
    Qhn(k) = abs(snm(k))*Tempr_p(j)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) - (0.5)*(I^2)*rn1(k); % - abs(0.5*(thnm(k))*I) ; % Heat input to hot side
    Qcn(k) = abs(snm(k))*Tempr_p(j+2)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) + 0.5*(I^2)*rn1(k); % + abs(0.5*((thnm(k)))*I) ;% Heat rejection from cold side
    Pn(k) = double(Qhn(k)- Qcn(k)); % Power output per leg


end
shp2 = hp2/hp2elems;
shn2 = hn2/hp2elems;

for k = 1:hp2elems
  
    if ( k ==1)
        j = C3*2 + 1;
    else
        j = C3*2 +1 + (k-1)*2;
    end
    
    
    
    knm(k) = (kn2(Tempr_p(j)) - kn2(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    kpm(k) = (kp2(Tempr_p(j)) - kp2(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));

    snm(k) =  (sn2(Tempr_n(j)) - sn2(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2));
    spm(k) =  (sp2(Tempr_p(j)) - sp2(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    ppm(k) =  (pp2(Tempr_p(j)) - pp2(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2));
    pnm(k) =  (pn2(Tempr_n(j)) - pn2(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2));
    thpm(k) =  (thp(Tempr_p(j)) - thp(Tempr_p(j+2)))/(Tempr_p(j)-Tempr_p(j+2)); % new code
    thpm(k) = double(thpm(k));
    thnm(k) =  (thn(Tempr_n(j)) - thn(Tempr_n(j+2)))/(Tempr_n(j)-Tempr_n(j+2)); % new code
    thnm(k) = double(thnm(k));

    rp12(k) = (ppm(k)* shp2 )/(lp*tp);  % [Ohms] Resistance p-leg
    rn12(k) = (pnm(k)* shn2 )/(ln*tn); % [Ohms] Resistance n-leg
    cp(k) =  (kpm(k)* lp*tp)/shp2; %  Thermal Conductance of p-leg
    cn(k) = (knm(k)* ln*tn)/shn2; % Thermal Conductance of n-leg
    
    V_oc2_e(k) = (spm(k)-snm(k))*(Tempr_p(j) -Tempr_p(j+2));
    
    Qhp2(k) = spm(k)*Tempr_p(j)*I  + cp(k)*(Tempr_p(j)-Tempr_p(j+2))-(0.5)*(I^2)*rp12(k) ; % - abs(0.5*thpm(k)*I);  % Heat input to hot side
    Qcp2(k) = spm(k)*Tempr_p(j+2)*I + cp(k)*(Tempr_p(j)-Tempr_p(j+2))+ 0.5*(I^2)*rp12(k) ; %+ abs(0.5*thpm(k)*I); % Heat rejection from cold side
    Pp2(k) = double((Qhp2(k)- Qcp2(k))) ;% Power output per leg
   
    Qhn2(k) = abs(snm(k))*Tempr_p(j)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) - (0.5)*(I^2)*rn12(k); % - abs(0.5*(thnm(k))*I) ; % Heat input to hot side
    Qcn2(k) = abs(snm(k))*Tempr_p(j+2)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) + 0.5*(I^2)*rn12(k); % + abs(0.5*((thnm(k)))*I) ;% Heat rejection from cold side
    Pn2(k) = double(Qhn2(k)- Qcn2(k)); % Power output per leg


end

% Qhin1p = max(Qhp)
% Qhin1n =  max(Qhn);
% Qhin2p = max(Qhp);
% Qhin2n = max(Qhn);
Qhin = max([Qhp, Qhp2]) + max([Qhn, Qhn2])
Rpl = sum(rp1);
Rnl = sum(rn1);
Rpl2 = sum(rp12);
Rnl2 = sum(rn12);

R_FE = Rpl + Rnl + Rpl2 + Rnl2 +  Rt1t + Rb1tn + Rb1tp +CR;
TotalP = sum(Pp);
TotalN = sum(Pn);
TotalP2 = sum(Pp2);
TotalN2 = sum(Pn2);
Power = TotalP + TotalN + TotalP2 + TotalN2 - (Rt1t + Rb1tn + Rb1tp + CR)*I^2
V_oc = sum(V_oc_e) + sum(V_oc2_e)
I_voc = V_oc/(R_FE*2)
Power_El = (V_oc^2)/(4*R_FE)
Eff = Power_El*100/Qhin
RLoad = linspace (0, 3.5*R_FE, 1000);


P_swipe = ((V_oc.^2).*(RLoad))./((RLoad+R_FE).^2);
P_el_max = max(P_swipe);
Eff_swipe = 100*P_swipe./Qhin;
[ max_Eff, Ind] = max(Eff_swipe)
Pow_maxeff = P_swipe(Ind)
RLoad_maxeff = RLoad(Ind)
I_swipe = V_oc./(R_FE+RLoad);

VTEG = R_FE.*(V_oc)./(RLoad + R_FE);
ITEG = V_oc./(RLoad + R_FE);
ITEG(end +1) = 0;
VTEG = ITEG*R_FE;
Vdev = V_oc - (R_FE*ITEG);
T_LHot = Tempr_p(C2*2 +1)
T_LInt = Tempr_p(C3*2 +1)
T_LCold = Tempr_p(C4*2 +1)
toc


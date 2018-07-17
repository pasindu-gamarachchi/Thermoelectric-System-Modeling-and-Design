%% Finite Element Model for TE Unicouple
% Pasindu Gamarachchi - Email: pgamarachchi@gmail.com
% Input hot and cold side temperatures, along with unicouple dimensions and
% temperature dependent material properties.

clc 
clear all

% Naming Convention
% property/leg/otherconsideration   no slashes

% Nodes and Elements
elems =100; % 
numn = 2*elems +1;

% Input Temperatures
Th = 600;
Tc = 100;   


% Input Dimensions

% N - Leg
hn = (2.4*10^-3);
ln = 1.5*10^-3; 
tn = 1.5*10^-3;
% P - Leg

hp = (2.4*10^-3); 
lp = 1.5*10^-3; 
tp = 1.5*10^-3 ; 


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

% Input Temperature Dependent properties
syms T


% Seebeck - N 
TNS = [ 29.85	49.85	99.85	149.85	199.85	249.85	299.85	349.85	399.85	449.85	499.85	549.85	599.85	649.85	651];
SN = [ -0.000315925	-0.000318116	-0.000318488	-0.000312653	-0.000301963	-0.000287772	-0.000271431	-0.000254295	-0.000237715	-0.000223045	-0.000211637	-0.000204845	-0.000204021	-0.000210518	-0.00021052];
CNS = polyfit( TNS, SN, 6);
snf = int ( CNS(1)*T^6 + CNS(2)*T^5 + CNS(3)*T^4 + CNS(4)*T^3 + CNS(5)*T^2 + CNS(6)*T + CNS(7));

% Electrical - N 
TNP = TNS;
PN = [ 0.000120408	0.000112708	9.49E-05	7.90E-05	6.51E-05	5.32E-05	4.33E-05	3.54E-05	2.95E-05	2.56E-05	2.37E-05	2.37E-05	2.58E-05	2.98E-05	2.98E-05];
CNP = polyfit( TNS, PN, 6);
pnf = int ( CNP(1)*T^6 + CNP(2)*T^5 + CNP(3)*T^4 + CNP(4)*T^3 + CNP(5)*T^2 + CNP(6)*T + CNP(7));

% Theramal - N 
TNK = TNS;
KN = [1.69829	1.62797	1.47154	1.34126	1.23517	1.15116	1.08704	1.04049	1.00908	0.99027	0.9814	0.97972	0.98233	0.98626	0.98626];
CNK = polyfit( TNK, KN, 6);
knf = int ( CNK(1)*T^6 + CNK(2)*T^5 + CNK(3)*T^4 + CNK(4)*T^3 + CNK(5)*T^2 + CNK(6)*T + CNK(7));

% Seebeck - P 
TPS = [ 27.3	49.1	99.3	150	200.1	250.1	299.8	350	400.2	450.1	500.2	550.1	600 	650];
SP = [9.09E-05	0.000100049	0.000127219	0.000154621	0.000182598	0.000207833	0.000234685	0.000251779	0.000269364	0.000278921	0.000288352	0.000292903	0.000296377	0.000295786 ];
CPS = polyfit( TPS, SP, 6);
spf = int ( CPS(1)*T^6 + CPS(2)*T^5 + CPS(3)*T^4 + CPS(4)*T^3 + CPS(5)*T^2 + CPS(6)*T + CPS(7));

% Electrical - P
TPP = TPS;
PP = [4.90E-06	5.30E-06	6.53E-06	8.38E-06	1.08E-05	1.39E-05	1.74E-05	2.02E-05	2.28E-05	2.47E-05	2.73E-05	2.97E-05	3.13E-05	3.26E-05 ];
CPP = polyfit( TPP, PP, 6);
ppf = int ( CPP(1)*T^6 + CPP(2)*T^5 + CPP(3)*T^4 + CPP(4)*T^3 + CPP(5)*T^2 + CPP(6)*T + CPP(7));

% Thermal - P
TPK = TPS;
KP = [2.9113	2.74663	2.37698	2.03981	1.76258	1.53701	1.33858	1.22853	1.13914	1.09221	1.04702	1.02015	1.01201	1.01397];
CPK = polyfit( TPK, KP, 6);
kpf = int ( CPK(1)*T^6 + CPK(2)*T^5 + CPK(3)*T^4 + CPK(4)*T^3 + CPK(5)*T^2 + CPK(6)*T + CPK(7));
 
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


kn = symfun( knf, T); % Temperature dependent k for n-leg
kp = symfun (kpf, T); % Temperature dependent k for p-leg
sn = symfun (snf, T); % Temperature dependent s for n-leg
sp = symfun (spf, T); % Temperature dependent k for p-leg
pp = symfun (ppf, T); % Temperature dependent k for p-leg
pn = symfun (pnf, T); % Temperature dependent k for p-leg
kc = symfun (kcf, T); % Temperature dependent k for Copper 102
pc = symfun (pcf, T); % Temperature dependent p for Copper 102
ka = symfun (kaf, T); % Temperature dependent k for Al2O3


% Integral Averages for initial Temperature profile

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

kcm =  (kc(Th) - kc(Tc))/(Th-Tc);
kcm = double(kcm);
pcm =  (pc(Th) - pc(Tc))/(Th-Tc);
pcm = double(pcm);

kam =  (ka(Th) - ka(Tc))/(Th-Tc);
kam = double(kam);




% Thermal Circuit for initial temperature profile
rc1 = (t3t/(kcm*t3l*t3w/2));
ra1 = (t2t/(kam*t2l*t2w/2));
rc2 = (t1t/(kcm*t1l*t1w/2));
cp =  (kpm* lp*tp)/hp;  
rthp = 1/cp;
cn =  (knm* ln*tn)/hn; 
rthn = 1/cn;
RLegs = 1/(cp +cn);

ra2 = (b2t/(kam.*b2l*b2w/2));
rc3 = (b1t/(kcm*b1l*b1w/2));
rc4 = (b3t/(kcm*b3l*b3w/2));

Rtot = rc1 + ra1 + rc2 + RLegs + rc3 + ra2 + rc4;

P1 = rc1/Rtot;
P2 = ra1/Rtot;
P3 = rc2/Rtot;
P4 = RLegs/Rtot;
P5 = rc3/Rtot;
P6 = ra2/Rtot;
P7 = rc4/Rtot;

DelTP1 =  (P1*(Th-Tc));
DelTP2 = (P2*(Th-Tc));
DelTP3 = (P3*(Th-Tc));
DelTP4 = (P4*(Th-Tc));
DelTP5 = (P5*(Th-Tc));
DelTP6 = (P6*(Th-Tc));
DelTP7 = (P7*(Th-Tc));


DeltT = DelTP1+ DelTP2 + DelTP3+ DelTP4+ DelTP5+ DelTP6 + DelTP7;

Tcu1 = Th - DelTP1;
Ta1 = Tcu1 - DelTP2;
Tcu2 = Ta1 - DelTP3;
Tpleg = Tcu2 - DelTP4;
Tcu3 = Tpleg - DelTP5;
Ta2 = Tcu3 - DelTP6;
Tcu4 = Ta2 - DelTP7;

Tp = [ Th Tcu1 Ta1 Tcu2 Tpleg Tcu3 Ta2  Tc];

%%
TotalHeight = t3t + hp + t2t + t1t + b3t + b2t + b1t;

t3telems = round((t3t/TotalHeight)*elems);
t2telems = round((t2t/TotalHeight)*elems);
t1telems = round((t1t/TotalHeight)*elems);
hpelems = round((hp/TotalHeight)*elems);
b3telems = round((b3t/TotalHeight)*elems);
b2telems = round((b2t/TotalHeight)*elems);
b1telems = round((b1t/TotalHeight)*elems);
summedelems = sum([t3telems , t2telems, t1telems, hpelems, b3telems, b2telems ,b1telems]);

while ( summedelems ~= elems)

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


for i = 1: t3telems 
    
    dt = (Th - Tcu1)/t3telems;
    T_t3t(i+1) = Th - dt*(i);
    
end
T_t3t(1) = Th;

for i = 1: t2telems 
    
    dt = (Tcu1 - Ta1)/t2telems;
    T_t2t(i+1) = Tcu1 - dt*(i);
    
end
T_t2t(1) = Tcu1;

for i = 1: t1telems 
    
    dt = (Ta1 - Tcu2)/t1telems;
    T_t1t(i+1) = Tcu2 - dt*(i);
    
end
T_t1t(1) = Ta1;

for i = 1: hpelems 
    
    dt = (Tcu2 - Tpleg)/hpelems;
    T_hp(i+1) = Tcu2 - dt*(i);
    
end
T_hp(1) = Tcu2;

for i = 1: b1telems 
    
    dt = (Tpleg - Tcu3)/b1telems;
    T_b1t(i+1) = Tpleg - dt*(i);
    
end
T_b1t(1) = Tpleg;


for i = 1: b2telems 
    
    dt = (Tcu3 - Ta2)/b2telems;
    T_b2t(i+1) = Tcu3 - dt*(i);
    
end
T_b2t(1) = Tcu3;

for i = 1: b3telems 
    
    dt = (Ta2 - Tc)/b3telems;
    T_b3t(i+1) = Ta2 - dt*(i);
    
end
T_b3t(1) = Ta2;
T_b3t(end) = Tc;


Tp = [ T_t3t T_t2t T_t1t T_hp T_b1t T_b2t T_b3t ];
%%

ppm =  (pp(T_hp(1)) - pp(T_hp(end)))/(T_hp(1)-T_hp(end));
ppm = double(ppm);
pnm =  (pn(T_hp(1)) - pn(T_hp(end)))/(T_hp(1)-T_hp(end));
pnm = double(pnm);

rp = (ppm * hp )/(lp*tp); 
rn = (pnm * hn )/(ln*tn); 

pt1tm =  (pc(T_t1t(1)) - pc(T_t1t(end)))/(T_t1t(1) - T_t1t(end) );
pt1tm = double(pt1tm);

pb1tpm =  (pc(T_b1t(1)) - pc(T_b1t(end)))/(T_b1t(1) - T_b1t(end) );
pb1tpm = double(pb1tpm);

pb1tnm =  (pc(T_b1t(1)) - pc(T_b1t(end)))/(T_b1t(1) - T_b1t(end) );
pb1tnm = double(pb1tnm);



snm =  (sn(T_hp(1)) - sn(T_hp(end)))/(T_hp(1)-T_hp(end));
snm = double(snm);
spm =  (sp(T_hp(1)) - sp(T_hp(end)))/(T_hp(1)-T_hp(end));
spm = double(spm);

% Contact Resistance
C_pho = 10*(10^-6);
CR_n = C_pho/(tp*100)*(lp*100);
CR_p = C_pho/(tp*100)*(lp*100);
CR = 2*(CR_n + CR_p);


Rt1t = (pt1tm* t1w)/(t1l*t1t);
Rb1tp = (pb1tpm* b1w)/(b1l*b1t);
Rb1tn =  (pb1tnm* b1w)/(b1l*b1t) ;
Rt = rp +rn + Rt1t + Rb1tp +Rb1tn + CR;
RL = Rt ; 
R = RL + Rt ;  
S =  spm -snm ; 
I = (S*(T_hp(1) -T_hp(end)))/(R);

%%
t1te_h = t1t/t1telems;

Jh_t1t = (I^2)* Rt1t;
Jh_b1tp = (I^2)*Rb1tp;
Jh_b1tn = (I^2)*Rb1tn;

Jh_t1t_e = Jh_t1t/(t1telems*2);
Jh_b1tp_e = Jh_t1t/(b1telems);
Jh_b1tn_e = Jh_t1t/(b1telems);

shp = hp/hpelems;
shn = hn/hpelems;

%% Thermoelectric Energy Generation calculations using initial temperature profile  guess 
%% to be input into FE calculations
for j = 1:hpelems
  
    knm(j) = (kn(T_hp(j)) - kn(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));
    kpm(j) = (kp(T_hp(j)) - kp(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));

    snm(j) =  (sn(T_hp(j)) - sn(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));
    spm(j) =  (sp(T_hp(j)) - sp(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));
    ppm(j) =  (pp(T_hp(j)) - pp(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));
    pnm(j) =  (pn(T_hp(j)) - pn(T_hp(j+1)))/(T_hp(j)-T_hp(j+1));
  

    rp1(j) = (ppm(j)* shp )/(lp*tp); 
    rn1(j) = (pnm(j)* shn )/(ln*tn); 
    cp(j) =  (kpm(j)* lp*tp)/shp; 
    cn(j) = (knm(j)* ln*tn)/shn; 
    
    Qhp(j) = spm(j)*T_hp(j)*I + cp(j)*(T_hp(j)-T_hp(j+1))-(0.5)*(I^2)*rp1(j); 
    Qcp(j) = spm(j)*T_hp(j+1)*I + cp(j)*(T_hp(j)-T_hp(j+1))+ 0.5*(I^2)*rp1(j); 
    Pp(j) = double((Qhp(j)- Qcp(j))) ;
   
    Qhn(j) = abs(snm(j))*T_hp(j)*I + cn(j)*(T_hp(j)-T_hp(j+1)) - (0.5)*(I^2)*rn1(j) ; 
    Qcn(j) = abs(snm(j))*T_hp(j+1)*I + cn(j)*(T_hp(j)-T_hp(j+1)) + 0.5*(I^2)*rn1(j) ; 
    Pn(j) = double(Qhn(j)- Qcn(j)); 


end


% Temperature Dependent elemental k - matrix
for i = 1:t3telems 
    
    k_t3tp_e(i) = double((kc(T_t3t(i)) - kc(T_t3t(i+1)))/(T_t3t(i) - T_t3t(i+1)));
    k_t3tn_e(i) = double((kc(T_t3t(i)) - kc(T_t3t(i+1)))/(T_t3t(i) - T_t3t(i+1)));
    A_t3t(i) = t3w*t3l;

end

for i = 1:t2telems 
    
    k_t2tp_e(i) = double((ka(T_t2t(i)) - ka(T_t2t(i+1)))/(T_t2t(i) - T_t2t(i+1)));
    k_t2tn_e(i) = double((ka(T_t2t(i)) - ka(T_t2t(i+1)))/(T_t2t(i) - T_t2t(i+1)));
    A_t2tp(i) = 0.5*(t2w*t2l);
    A_t2tn(i) = 0.5*(t2w*t2l);

end

for i = 1:t1telems 
    
    k_t1tp_e(i) = double((kc(T_t1t(i)) - kc(T_t1t(i+1)))/(T_t1t(i) - T_t1t(i+1)));
    k_t1tn_e(i) = double((kc(T_t1t(i)) - kc(T_t1t(i+1)))/(T_t1t(i) - T_t1t(i+1)));
    A_t1tp(i) = 0.5*(t1w*t1l);
    A_t1tn(i) = 0.5*(t1w*t1l);
end

for i = 1:hpelems 
    
    k_hp_e(i) = double((kp(T_hp(i)) - kp(T_hp(i+1)))/(T_hp(i) - T_hp(i+1)));
    k_hn_e(i) = double((kn(T_hp(i)) - kn(T_hp(i+1)))/(T_hp(i) - T_hp(i+1)));
    A_p(i) = lp*tp;
    A_n(i) = ln*tn;
    
    
end

for i = 1:b1telems 
    
    k_b1tp_e(i) = double((kc(T_b1t(i)) - kc(T_b1t(i+1)))/(T_b1t(i) - T_b1t(i+1)));
    k_b1tn_e(i) = double((kc(T_b1t(i)) - kc(T_b1t(i+1)))/(T_b1t(i) - T_b1t(i+1)));
    A_b1t(i) = b1w*b1l;
end

for i = 1:b2telems 
    
    k_b2tp_e(i) = double((ka(T_b2t(i)) - ka(T_b2t(i+1)))/(T_b2t(i) - T_b2t(i+1)));
    k_b2tn_e(i) = double((ka(T_b2t(i)) - ka(T_b2t(i+1)))/(T_b2t(i) - T_b2t(i+1)));
    A_b2tp(i) = 0.5*(b2w*b2l);
    A_b2tn(i) = 0.5*(b2w*b2l);
end


for i = 1:b3telems 
    
    k_b3tp_e(i) = double((kc(T_b3t(i)) - kc(T_b3t(i+1)))/(T_b3t(i) - T_b3t(i+1)));
    k_b3tn_e(i) = double((kc(T_b3t(i)) - kc(T_b3t(i+1)))/(T_b3t(i) - T_b3t(i+1)));
    A_b3tp(i) = 0.5*(b3w*b3l);
    A_b3tn(i) = 0.5*(b3w*b3l);

end

kp_elems = [ k_t3tp_e, k_t2tp_e, k_t1tp_e, k_hp_e, k_b1tp_e, k_b2tp_e, k_b3tp_e];
kn_elems = [ k_t3tn_e, k_t2tn_e, k_t1tn_e, k_hn_e, k_b1tn_e, k_b2tn_e, k_b3tn_e];

Ap_elems = [ A_t3t, A_t2tp, A_t1tp, A_p, A_b1t, A_b2tp, A_b3tp];
An_elems = [ A_t3t, A_t2tn, A_t1tn, A_n, A_b1t, A_b2tn, A_b3tn];

Kp = sparse(numn, numn);
Kn = sparse(numn, numn);

el_l = TotalHeight/elems;

%% Global K- Matrix Assembly
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

% F - Vector assembly using TE energy generation and joule heating
C0 = t3telems;
C1 = t3telems + t2telems; 
C2 = C1 + t1telems; 
C3 = C2 + hpelems; 
C4 = C3 + b1telems; 
C5 = C4 + b2telems;

for i = C1+1: C2
    
    E_gen_t1t(i) = Jh_t1t_e/(Ap_elems(i) * el_l);
    
end

for i = 1:hpelems

    E_genp(i) = -Pp(i)/(Ap_elems(i+C2)*el_l);
    E_genn(i) = -Pn(i)/(Ap_elems(i+C2)*el_l);

end

for i = C3+1: C4
    
    j= i -C3;
    E_gen_b1tp(j) = Jh_b1tp_e/(Ap_elems(i)* el_l);
    E_gen_b1tn(j) = Jh_b1tn_e/(Ap_elems(i)* el_l);

end
E_gen_b2t = zeros(b2telems,1)';
E_gen_b3t = zeros(b3telems,1)';

E_Gen_pL = [ E_gen_t1t, E_genp, E_gen_b1tp, E_gen_b2t, E_gen_b3t];
E_Gen_nL = [ E_gen_t1t, E_genn, E_gen_b1tn, E_gen_b2t, E_gen_b3t];
    
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

%% Adjusing for Boundary Conditions
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

F_nold = F_n;
F_pold = F_p;


ppm =  (pp(Tempr_p(C2*2 +1)) - pp(Tempr_p(C3*2 +1)))/(Tempr_p(C2*2 +1)-Tempr_p(C3*2 +1));
ppm = double(ppm);
pnm =  (pn(Tempr_p(C2*2 +1)) - pn(Tempr_p(C3*2 +1)))/(Tempr_p(C2*2 +1)-Tempr_p(C3*2 +1));
pnm = double(pnm);

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

Rt1t = (pt1tm* t1w)/(t1l*t1t);
Rb1tp = (pb1tpm* b1w)/(b1l*b1t);
Rb1tn =  (pb1tnm* b1w)/(b1l*b1t) ;
Rt = rp +rn + Rt1t + Rb1tp +Rb1tn + CR ;
RL = Rt ; 
R = RL + Rt ; 

S =  spm -snm ;
Voc_av = S*(Tempr_p(C2*2 +1) -Tempr_p(C3*2 +1));
I = S*(Tempr_p(C2*2 +1) -Tempr_p(C3*2 +1))./(R);


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

    rp1(k) = (ppm(k)* shp )/(lp*tp); 
    rn1(k) = (pnm(k)* shn )/(ln*tn); 
    cp(k) =  (kpm(k)* lp*tp)/shp; 
    cn(k) = (knm(k)* ln*tn)/shn; 
    
    V_oc_e(k) = (spm(k)-snm(k))*(Tempr_p(j) -Tempr_p(j+2));
    
    Qhp(k) = spm(k)*Tempr_p(j)*I  + cp(k)*(Tempr_p(j)-Tempr_p(j+2))-(0.5)*(I^2)*rp1(k) ; 
    Qcp(k) = spm(k)*Tempr_p(j+2)*I + cp(k)*(Tempr_p(j)-Tempr_p(j+2))+ 0.5*(I^2)*rp1(k) ; 
    Pp(k) = double((Qhp(k)- Qcp(k))) ;
   
    Qhn(k) = abs(snm(k))*Tempr_p(j)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) - (0.5)*(I^2)*rn1(k); 
    Qcn(k) = abs(snm(k))*Tempr_p(j+2)*I + cn(k)*(Tempr_p(j)-Tempr_p(j+2)) + 0.5*(I^2)*rn1(k); 
    Pn(k) = double(Qhn(k)- Qcn(k)); 

end

%% Iterative Process for Temperature profiles to converge
it =1;
cc = 1;
err = cc*3; % Initalize Error value
while err > cc

    if it ==1 
        Tempr_p = Tempr_p;
        Tempr_n = Tempr_n;
    else
        Tempr_p = Tempr_p2;
        Tempr_n = Tempr_n2;
    end
%% Temperature Dependent elemental k - matrix using new Temperature profile
for i = 1:t3telems %  t3t (Very Top Cu Layer) Copper Conductivity and Area
    
    if ( i ==1)
        j =  1;
    else
        j = (i-1)*2 +1;
    end
    
    k_t3tp_e(i) = double((kc(Tempr_p(j)) - kc(Tempr_p(j+2)))/(Tempr_p(j) - Tempr_p(j+1)));
    k_t3tn_e(i) = double((kc(Tempr_n(j)) - kc(Tempr_n(j+2)))/(Tempr_n(j) - Tempr_n(j+2)));
  
end

for i = 1:t2telems % t22 Top Alumina Conductivity and Area
    
     if ( i ==1)
        j = C0 + 1;
    else
        j = C0 +(i-1)*2 +1;
    end
    
    
     k_t2tp_e(i) = double((ka(Tempr_p(i)) - ka(Tempr_p(i+1)))/(Tempr_p(i) - Tempr_p(i+1)));
     k_t2tn_e(i) = double((ka(Tempr_n(i)) - ka(Tempr_n(i+1)))/(Tempr_n(i) - Tempr_n(i+1)));


end

for i = 1:t1telems % t1t ( Copper Layer Connecting Legs ) Conductivity and Area
    
     if ( i ==1)
        j = C1 + 1;
    else
        j = C1 +(i-1)*2 +1;
    end
    
     k_t1tp_e(i) = double((kc(Tempr_p(i)) - kc(Tempr_p(i+1)))/(Tempr_p(i) - Tempr_p(i+1)));
     k_t1tn_e(i) = double((kc(Tempr_n(i)) - kc(Tempr_n(i+1)))/(Tempr_n(i) - Tempr_n(i+1)));
  
end

for i = 1:hpelems % Leg conductivity and area
    
     if ( i ==1)
        j = C2 + 1;
    else
        j = C2 +(i-1)*2 +1;
    end
    
    k_hp_e(i) = double((kp(Tempr_p(i)) - kp(Tempr_p(i+1)))/(Tempr_p(i) - Tempr_p(i+1)));
     k_hn_e(i) = double((kn(Tempr_n(i)) - kn(Tempr_n(i+1)))/(Tempr_n(i) - Tempr_n(i+1)));

    
end

for i = 1:b1telems % Bottom Copper Conductivity and Area
    
    if ( i ==1)
        j = C3 + 1;
    else
        j = C3 +(i-1)*2 +1;
    end
    
     k_b1tp_e(i) = double((kc(Tempr_p(i)) - kc(Tempr_p(i+1)))/(Tempr_p(i) - Tempr_p(i+1)));
     k_b1tn_e(i) = double((kc(Tempr_n(i)) - kc(Tempr_n(i+1)))/(Tempr_n(i) - Tempr_n(i+1)));
    
end

for i = 1:b2telems % Bottom Alumina Conductivity and Area
    
    if ( i ==1)
        j = C3 + 1;
    else
        j = C3 +(i-1)*2 +1;
    end
    
    
     k_b2tp_e(i) = double((ka(Tempr_p(i)) - ka(Tempr_p(i+1)))/(Tempr_p(i) - Tempr_p(i+1)));
     k_b2tn_e(i) = double((ka(Tempr_n(i)) - ka(Tempr_n(i+1)))/(Tempr_n(i) - Tempr_n(i+1)));
 
end


for i = 1:b3telems % Bottomost Copper Conductivity and Area
    
    if ( i ==1)
        j = C4 + 1;
    else
        j = C4 +(i-1)*2 +1;
    end
    
     k_b3tp_e(i) = double((kc(Tempr_p(i)) - kc(Tempr_p(i+1)))/(Tempr_p(i) - Tempr_p(i+1)));
     k_b3tn_e(i) = double((kc(Tempr_n(i)) - kc(Tempr_n(i+1)))/(Tempr_n(i) - Tempr_n(i+1)));

end

kp_elems = [ k_t3tp_e, k_t2tp_e, k_t1tp_e, k_hp_e, k_b1tp_e, k_b2tp_e, k_b3tp_e];
kn_elems = [ k_t3tn_e, k_t2tn_e, k_t1tn_e, k_hn_e, k_b1tn_e, k_b2tn_e, k_b3tn_e];

Ap_elems = [ A_t3t, A_t2tp, A_t1tp, A_p, A_b1t, A_b2tp, A_b3tp];
An_elems = [ A_t3t, A_t2tn, A_t1tn, A_n, A_b1t, A_b2tn, A_b3tn];

Kp = sparse(numn, numn);
Kn = sparse(numn, numn);

el_l = TotalHeight/elems;

%% Global K- Matrix Assembly
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
F_p = sparse(1,numn);
F_n = sparse(1,numn);

% F - Vector assembly using TE energy generation and joule heating
C0 = t3telems;
C1 = t3telems + t2telems; 
C2 = C1 + t1telems; 
C3 = C2 + hpelems; 
C4 = C3 + b1telems; 
C5 = C4 + b2telems;

for i = C1+1: C2
    
  E_gen_t1t(i) = Jh_t1t_e/(Ap_elems(i) * el_l);
    
end

for i = 1:hpelems

    E_genp(i) = -Pp(i)/(Ap_elems(i+C2)*el_l);
    E_genn(i) = -Pn(i)/(Ap_elems(i+C2)*el_l);

end

for i = C3+1: C4
    
    j= i -C3;
    E_gen_b1tp(j) = Jh_b1tp_e/(Ap_elems(i)* el_l);
    E_gen_b1tn(j) = Jh_b1tn_e/(Ap_elems(i)* el_l);

end
E_gen_b2t = zeros(b2telems,1)';
E_gen_b3t = zeros(b3telems,1)';

E_Gen_pL = [ E_gen_t1t, E_genp, E_gen_b1tp, E_gen_b2t, E_gen_b3t];
E_Gen_nL = [ E_gen_t1t, E_genn, E_gen_b1tn, E_gen_b2t, E_gen_b3t];
    
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

%% Adjusing for Boundary Conditions
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


Tempr_p2 = Kp\F_p';
Tempr_n2 = Kn\F_n';

diffp = abs(Tempr_p -Tempr_p2) ;
diffn =  abs(Tempr_n -Tempr_n2);

err = sum(diffp + diffn);
it = it +1;

end


%% Thermoelectric Calculations using Final Temperature Profile
ppm =  (pp(Tempr_p2(C2*2 +1)) - pp(Tempr_p2(C3*2 +1)))/(Tempr_p2(C2*2 +1)-Tempr_p2(C3*2 +1));
ppm = double(ppm);
pnm =  (pn(Tempr_p2(C2*2 +1)) - pn(Tempr_p2(C3*2 +1)))/(Tempr_p2(C2*2 +1)-Tempr_p2(C3*2 +1));
pnm = double(pnm);

pt1tm =  (pc(Tempr_p2(C1*2 +1 )) - pc(Tempr_p2(C2*2 +1)))/(Tempr_p2(C1*2 +1) - Tempr_p2(C2*2 +1) );
pt1tm = double(pt1tm);

pb1tpm =  (pc(Tempr_p2(C3*2 +1)) - pc(Tempr_p2(C4*2 +1)))/(Tempr_p2(C3*2 +1) - Tempr_p2(C4*2 +1) );
pb1tpm = double(pb1tpm);

pb1tnm =  (pc(Tempr_n2(C3*2 +1)) - pc(Tempr_n2(C4*2 +1)))/(Tempr_n2(C3*2 +1) - Tempr_n2(C4*2 +1) );
pb1tnm = double(pb1tnm);

snm =  (sn(Tempr_n2(C2*2 +1)) - sn(Tempr_n2(C3*2 +1)))/(Tempr_n2(C2*2 +1)-Tempr_n2(C3*2 +1));
snm = double(snm);
spm =  (sp(Tempr_p2(C2*2 +1)) - sp(Tempr_p2(C3*2 +1)))/(Tempr_p2(C2*2 +1)-Tempr_p2(C3*2 +1));
spm = double(spm);

Rt1t = (pt1tm* t1w)/(t1l*t1t);
Rb1tp = (pb1tpm* b1w)/(b1l*b1t);
Rb1tn =  (pb1tnm* b1w)/(b1l*b1t) ;
Rt = rp +rn + Rt1t + Rb1tp +Rb1tn + CR ;
RL = Rt ;
R = RL + Rt ; 

S =  spm -snm ; 
Voc_av = S*(Tempr_p2(C2*2 +1) -Tempr_p2(C3*2 +1));
I = S*(Tempr_p2(C2*2 +1) -Tempr_p2(C3*2 +1))./(R);


for k = 1:hpelems
  
    if ( k ==1)
        j = C2*2 + 1;
    else
        j = C2*2 +1 + (k-1)*2;
    end
    
    
    
    knm(k) = (kn(Tempr_p2(j)) - kn(Tempr_p2(j+2)))/(Tempr_p2(j)-Tempr_p2(j+2));
    kpm(k) = (kp(Tempr_p2(j)) - kp(Tempr_p2(j+2)))/(Tempr_p2(j)-Tempr_p2(j+2));

    snm(k) =  (sn(Tempr_n2(j)) - sn(Tempr_n2(j+2)))/(Tempr_n2(j)-Tempr_n2(j+2));
    spm(k) =  (sp(Tempr_p2(j)) - sp(Tempr_p2(j+2)))/(Tempr_p2(j)-Tempr_p2(j+2));
    ppm(k) =  (pp(Tempr_p2(j)) - pp(Tempr_p2(j+2)))/(Tempr_p2(j)-Tempr_p2(j+2));
    pnm(k) =  (pn(Tempr_n2(j)) - pn(Tempr_n2(j+2)))/(Tempr_n2(j)-Tempr_n2(j+2));

    rp1(k) = (ppm(k)* shp )/(lp*tp);  
    rn1(k) = (pnm(k)* shn )/(ln*tn); 
    cp(k) =  (kpm(k)* lp*tp)/shp; 
    cn(k) = (knm(k)* ln*tn)/shn; 
    
    V_oc_e(k) = (spm(k)-snm(k))*(Tempr_p(j) -Tempr_p(j+2));
    
    Qhp(k) = spm(k)*Tempr_p2(j)*I  + cp(k)*(Tempr_p2(j)-Tempr_p2(j+2))-(0.5)*(I^2)*rp1(k) ; 
    Qcp(k) = spm(k)*Tempr_p2(j+2)*I + cp(k)*(Tempr_p2(j)-Tempr_p2(j+2))+ 0.5*(I^2)*rp1(k) ; 
    Pp(k) = double((Qhp(k)- Qcp(k))) ;
   
    Qhn(k) = abs(snm(k))*Tempr_p2(j)*I + cn(k)*(Tempr_p2(j)-Tempr_p2(j+2)) - (0.5)*(I^2)*rn1(k); 
    Qcn(k) = abs(snm(k))*Tempr_p2(j+2)*I + cn(k)*(Tempr_p2(j)-Tempr_p2(j+2)) + 0.5*(I^2)*rn1(k); 
    Pn(k) = double(Qhn(k)- Qcn(k)); 


end

Qhin = max(Qhp) + max(Qhn);
TotalP = sum(Pp);
TotalN = sum(Pn);
Power = TotalP + TotalN - (Rt1t + Rb1tn + Rb1tp + CR)*I^2;
V_oc = sum(V_oc_e);
Eff = Power*100/Qhin;


Qhin = max(Qhp) + max(Qhn);
Rpl = sum(rp1);
Rnl = sum(rn1);
R_FE = Rpl + Rnl + Rt1t + Rb1tn + Rb1tp +CR;
TotalP = sum(Pp);
TotalN = sum(Pn);
Power = TotalP + TotalN - (Rt1t + Rb1tn + Rb1tp + CR)*I^2
V_oc = sum(V_oc_e)
I_voc = V_oc/(R_FE*2)
Power_El = (V_oc^2)/(4*R_FE);
Eff = Power_El*100/Qhin
RLoad = linspace (0, 7.5*R_FE, 500);
% RLoad = linspace (0, 7.5*R_FE, 1000);


P_swipe = ((V_oc.^2).*(RLoad))./((RLoad+R_FE).^2);
P_el_max = max(P_swipe);
Eff_swipe = 100*P_swipe./Qhin;
[ max_Eff, Ind] = max(Eff_swipe)
Pow_maxeff = P_swipe(Ind);
RLoad_maxeff = RLoad(Ind);
I_swipe = V_oc./(R_FE+RLoad);

VTEG = R_FE.*(V_oc)./(RLoad + R_FE);
ITEG = V_oc./(RLoad + R_FE);
ITEG(end +1) = 0;
VTEG = ITEG*R_FE;
Vdev = V_oc - (R_FE*ITEG);
T_LHot = Tempr_p(C2*2 +1)
T_LCold = Tempr_p(C3*2 +1)

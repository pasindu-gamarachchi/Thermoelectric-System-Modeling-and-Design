%% Power Generation in a TE Couple, with Power output of 
% two legs seperated

% Updated 06/09
% DBC Dimensions : smaller dimensions from drawings
% Leg Properties : Half - Heusler measured 24/05 by Nick
% Leg Dimensions : Dimensions of unicouple made by Trent/Joe
% Updated 06/15
% Alumina and Copper Properties Updated

clc
clear 
close all
format long
tic

% Naming Convention
% property/leg/otherconsideration   no slashes

%  Th = [200 300 400 500 600 ];% [C] Hot side temperature
m = 400; % m is the number of segments, Set m >1
% RL = 0:0.001:0.2;
Thv = [ 200 300 400 500 600];

for p=1:length(Thv)


Th = Thv(p);
Tc = 100;   % [C] Cold side temperature
hn = 1.69*10^-3 ; % [m] n-leg height 
ln = 1.3*10^-3; %  [m] n-leg length
tn = 1.3*10^-3; %  [m] n-leg thickness
hp = 1.69*10^-3; %  [m] p-leg height
lp = 1.3*10^-3; % [m] p-leg length
tp = 1.3*10^-3 ; % [m] p-leg thickness
% NewI = 1.0;  %  Current Input, Load Resistance Set Manual
%  T3 - Copper
t3w = 1.93*10^-3;
t3l = 1.96*10^-3;
t3t = 0.2032*10^-3;

% T2 - Ceramic

t2w = 4.51*10^-3;
t2l = 2.26*10^-3;
t2t = 0.635*10^-3;

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

%% B1 - Copper

b1w = 4.07*10^-3;
b1l = 1.96*10^-3;
b1t = 0.2032*10^-3;

%% Temperature dependent properties/ Integrals of functions with respesct to T
% Electrical Resistivity N
syms x
pnf = int(3.79548552E-22*x^6 - 6.75630691E-19*x^5 + 4.62762403E-16*x^4 - 1.59187729E-13*x^3 + 2.48272530E-11*x^2 + 4.34512734E-09*x + 6.47172345E-06);

%% Seebeck N
snf = int(-3.07307824E-20*x^6 + 5.17273632E-17*x^5 - 3.14308248E-14*x^4 + 8.23370332E-12*x^3 - 6.95427569E-10*x^2 - 2.15931357E-07*x - 1.30728621E-04);

%%  Thermal Conductivity N
knf = int( -3.42296835E-16*x^6 + 5.14867564E-13*x^5 - 1.96473071E-10*x^4 - 4.33833033E-08*x^3 + 4.91278625E-05*x^2 - 1.38542448E-02*x + 5.56989213E+00);

%% Resisitivity P

ppf = int(1.26762914E-22*x^6 - 2.35719632E-19*x^5 + 1.60899624E-16*x^4 - 5.53161599E-14*x^3 + 1.53924136E-11*x^2 + 4.75622122E-09*x + 1.87165436E-06);


%% Seebeck P
 
 spf= int(-1.36579995E-20*x^6 + 2.24957017E-17*x^5 - 1.38487487E-14*x^4 + 3.96125625E-12*x^3 - 5.38600018E-10*x^2 + 2.14193134E-07*x + 6.98593256E-05);

 %% Thermal Conductivity P
 
 kpf = int( -8.25275365E-16*x^6 + 1.50160725E-12*x^5 - 1.02862881E-09*x^4 + 3.18113901E-07*x^3 - 3.14899647E-05*x^2 - 9.18620390E-03*x + 7.73547497E+00);
 

  %% Conductivity Copper 102 
 % New Code
 
 kcf = int(-2.82803550E-11*x^4 + 1.14240832E-07*x^3 - 1.55947483E-04*x^2 + 1.88171880E-02*x + 3.85115776E+02); % int( 3.98294713E-08*x^3 - 9.12017928E-05*x^2 - 2.16802231E-03*x + 3.87046692E+02);
 
 %% Resistivity Copper 102 
 % New Code
 pcf = int(1.59987141E-14*x^2 + 6.15360385E-11*x + 1.57380599E-08); 
 
 %% Conductivity Al2O3
 
 kaf = int(3.53561118E-11*x^4 - 1.04007564E-07*x^3 + 1.25189452E-04*x^2 - 7.77334115E-02*x + 2.62385634E+01); %int(1.79351954E-11*x^4 - 7.90633956E-08*x^3 + 1.30819330E-04*x^2 - 1.00219266E-01*x + 3.71764689E+01);
 
 %%
syms kn kp sn sp pp pn x % k: Conductivity, s: Seebeck, p: resistivity, p: p-leg, n:n-leg


kn = symfun( knf, x); % Temperature dependent k for n-leg
kp = symfun (kpf, x); % Temperature dependent k for p-leg
sn = symfun (snf, x); % Temperature dependent s for n-leg
sp = symfun (spf, x); % Temperature dependent k for p-leg
pp = symfun (ppf, x); % Temperature dependent k for p-leg
pn = symfun (pnf, x); % Temperature dependent k for p-leg
% New code
kc = symfun (kcf, x); % Temperature dependent k for Copper 102
pc = symfun (pcf, x); % Temperature dependent p for Copper 102
ka = symfun (kaf, x); % Temperature dependent k for Al2O3

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
% New code
kcm =  (kc(Th) - kc(Tc))./(Th-Tc);
kcm = double(kcm);
pcm =  (pc(Th) - pc(Tc))./(Th-Tc);
pcm = double(pcm);
%
kam =  (ka(Th) - ka(Tc))./(Th-Tc);
kam = double(kam);

rc1 = (t3t./(kcm.*t3l.*t3w));
ra1 = (t2t./(kam.*t2l.*t2w));
rc2 = (t1t./(kcm.*t1l.*t1w));
cp =  (kpm.* lp.*tp)./hp;  
cn =  (knm.* ln.*tn)./hn;  
rp = (ppm .* hp )./(lp.*tp);  % [Ohms] Resistance p-leg
rn = (pnm .* hn )./(ln.*tn); % [Ohms] Resistance n-leg



RC = (pcm.* t1w)./(t1l.*t1t);
Rt = rp +rn; % [Ohms] Resistance Total of legs
RBp = (pcm.* b1w)./(b1l.*b1t);
RBn = (pcm.* b1w)./(b1l.*b1t);
RL = Rt + RC +RBp +RBn; % [Ohms] Load resistance is set to internal resistance
R = RL +Rt +RC+RBp+RBn; % [Ohms] 
S =  spm -snm ; % Seebeck of couple
I = (S.*(Th -Tc))./(R);% [A] Current through the unicouple
% I = NewI;
Jhn = I.^2.*rn;
Jhp = I.^2.*rp;
%
%% Obtaining Tempr Profiles
Jhn = I.^2.*rn;
Jhp = I.^2.*rp;
Jhc = I.^2.* RC;
Jhbn = I.^2.*RBp;
Jhbp = I.^2.*RBn;

H = t1t+ t2t+t3t+b1t+b2t+b3t+hp;

% Thermal Resistances
rc3 = (b1t./(kcm.*b1l.*b1w));
ra2 = (b2t./(kam.*b2l.*b2w));
rc4 = (b3t./(kcm.*b3l.*b3w));
rtp = (1./cp);
Rtot = (rc1 + ra1 +rc2+(rtp)+rc3 +ra2+rc4);
A1 = (rc1./(Rtot));
A2 = ra1./(Rtot);
A3 = rc2./Rtot;
A4 = rtp./(Rtot);
A5 = rc3./Rtot;
A6 = ra2./Rtot;
A7 = rc4./Rtot;

DelT1 =  (A1.*(Th-Tc));
DelT2 = (A2.*(Th-Tc));
DelT3 = (A3.*(Th-Tc));
DelT4 = (A4.*(Th-Tc));
DelT5 = (A5.*(Th-Tc));
DelT6 = (A6.*(Th-Tc));
DelT7 = (A7.*(Th-Tc));

DeltT = DelT1+ DelT2 + DelT3+ DelT4+ DelT5+ DelT6 + DelT7;

Tcu1 = Th - DelT1;
Tce1 = Tcu1 - DelT2;
Tcu2 = Tce1 - DelT3;
Tpleg = Tcu2 - DelT4;
Tcu3 = Tpleg - DelT5;
Tce2 = Tcu3 - DelT6;
Tcu4 = Tce2 - DelT7;

T = [ Th Tcu1 Tce1 Tcu2 Tpleg Tcu3 Tce2  Tc];
% y = 0.0097978;
% z = [ y (y-t3t) (y-t3t-t2t)  (y-t3t-t2t-t1t) (y-t3t-t2t-t1t-hp) (y-t3t-t2t-t1t-hp-b3t) (y-t3t-t2t-t1t-hp-b3t-b2t) (y-t3t-t2t-t1t-hp-b3t-b2t-b1t) ];

m1 = 0.05*m;
m2 = 0.15*m;
m3 = 0.05*m;
m4 = 0.5*m;
m5 = 0.05*m;
m6 = 0.15*m;
m7 = 0.05*m;
sht3 = t3t./m1;

e =ones(m1,1);
AT1 = (1/sht3.^2).*spdiags( [e -2.*e e], -1:1, m1-1, m1-1);
FT1 = zeros(m1-1,1);
FT1(1) = -T(1)./(sht3.^2);
FT1(end)= -T(2)./(sht3.^2);

TT1 = AT1\FT1;
TT1(2:end+1)= TT1;
TT1(1) = T(1);

%%
sht2 = t2t./m2;
e =ones(m2+1,1);
AT2 = (1/sht2.^2).*spdiags( [e -2.*e e], -1:1, m2-1, m2-1);
FT2 = zeros(m2-1,1);
FT2(1) = -T(2)./(sht2.^2);
FT2(end)= -T(3)./(sht2.^2);

TT2 = AT2\FT2;
TT2(2:end+1)= TT2;
TT2(1) = T(2);

sht1 = t1t./m3;
e =ones(m3+1,1);
AT3 = (1/sht1.^2).*spdiags( [e -2.*e e], -1:1, m3-1, m3-1);
FT3 = (-Jhc./kcm).*ones(m3-1,1) ;
FT3(1) = -T(3)./(sht1.^2);
FT3(end)= -T(4)./(sht1.^2);

TT3 = AT3\FT3;
TT3(2:end+1)= TT3;
TT3(1) = T(3);


sht4 = hp./m4;
e =ones(m4+1,1);
AT4 = (1/sht4.^2).*spdiags( [e -2.*e e], -1:1, m4-1, m4-1);
FT4 = (-Jhp./kpm).*ones(m4-1,1) ;
FT4(1) = -T(4)./(sht4.^2);
FT4(end)= -T(5)./(sht4.^2);

TT4 = AT4\FT4;
TT4(2:end+1)= TT4;
TT4(1) = T(4);

sht5 = b1t./m5;
e =ones(m5+1,1);
AT5 = (1/sht5.^2).*spdiags( [e -2.*e e], -1:1, m5-1, m5-1);
FT5 = (-Jhc./kcm).*ones(m5-1,1) ;
FT5(1) = -T(5)./(sht5.^2);
FT5(end)= -T(6)./(sht5.^2);

TT5 = AT5\FT5;
TT5(2:end+1)= TT5;
TT5(1) = T(5);

sht6 = b2t./m6;
e =ones(m6+1,1);
AT6 = (1/sht6.^2).*spdiags( [e -2.*e e], -1:1, m6-1, m6-1);
FT6 = zeros(m6-1,1) ;
FT6(1) = -T(6)./(sht6.^2);
FT6(end)= -T(7)./(sht6.^2);

TT6 = AT6\FT6;
TT6(2:end+1)= TT6;
TT6(1) = T(6);

sht6 = b2t./m6;
e =ones(m6+1,1);
AT6 = (1/sht6.^2).*spdiags( [e -2.*e e], -1:1, m6-1, m6-1);
FT6 = zeros(m6-1,1) ;
FT6(1) = -T(6)./(sht6.^2);
FT6(end)= -T(7)./(sht6.^2);

TT6 = AT6\FT6;
TT6(2:end+1)= TT6;
TT6(1) = T(6);

sht7 = b2t./m7;
e =ones(m7+1,1);
AT7 = (1/sht7.^2).*spdiags( [e -2.*e e], -1:1, m7-1, m7-1);
FT7 = zeros(m7-1,1) ;
FT7(1) = -T(7)./(sht7.^2);
FT7(end)= -T(8)./(sht7.^2);

TT7 = AT7\FT7;
TT7(2:end+1)= TT7;
TT7(1) = T(7);
TT7(end+1) = Tc;

ATemp = [ TT1; TT2; TT3; TT4; TT5; TT6; TT7];


 for i = 1:m1 % T3 Copper
     
    kt1(i) = (kc(ATemp(i)) - kc(ATemp(i+1)))./((ATemp(i) -ATemp(i+1)));
    
 end
kt1 = double(kt1);
for i = 1:m2 % T2 Al2O3
    
    kt2(i) = (ka(ATemp(i+m1)) - ka(ATemp(i+m1+1)))./((ATemp(i+m1) -ATemp(i+m1+1)));
    
end
 kt2 = double(kt2);
for i = 1:m3 % T1 Copper
    
    kt3(i) = (kc(ATemp(i+(m1+m2))) - kc(ATemp(i+(m1+m2)+1)))./((ATemp(i+(m1+m2)) -ATemp(i+(m1+m2)+1)));
    
end
kt3 = double(kt3);
for i = 1:m4 % T1 Copper
    
    kt4(i) = (kp(ATemp(i+(m1+m2+m3))) - kp(ATemp(i+(m1+m2+m3)+1)))./((ATemp(i+(m1+m2+m3)) -ATemp(i+(m1+m2+m3)+1)));
    
end
kt4 = double(kt4);
for i = 1:m5 % T1 Copper
    
    kt5(i) = (kc(ATemp(i+(m1+m2+m3+m4))) - kc(ATemp(i+(m1+m2+m3+m4)+1)))./((ATemp(i+(m1+m2+m3+m4)) -ATemp(i+(m1+m2+m3+m4)+1)));
    
end
kt5 = double(kt5);
for i = 1:m6 % T1 Copper
    
    kt6(i) = (ka(ATemp(i+(m1+m2+m3+m4+m5))) - ka(ATemp(i+(m1+m2+m3+m4+m5)+1)))./((ATemp(i+(m1+m2+m3+m4+m5)) -ATemp(i+(m1+m2+m3+m4+m5)+1)));
    
end
kt6 = double(kt6);
for i = 1:m7 % T1 Copper
    
    kt7(i) = (kc(ATemp(i+(m1+m2+m3+m4+m5+m6))) - kc(ATemp(i+(m1+m2+m3+m4+m5+m6)+1)))./((ATemp(i+(m1+m2+m3+m4+m5+m6)) -ATemp(i+(m1+m2+m3+m4+m5+m6)+1)));
    
end
kt7 = double(kt7);

kt = [kt1 kt2 kt3 kt4 kt5 kt6 kt7]';

% Mid-point K's

for i=2:m-1
    kh(i) = (kt(i+1)+ kt(i))./2;
 
end

kh(1) = kt(1);
kh(end+1) = kt(end);
kh(end+1) = kt(end);


for i =1:m-1
    d1(i) = kh(i+1);
end

for i =1:m
    d2(i) = -(kh(i)+ kh(i+1));
end
C0 = 0.05*m;
C1 = C0 +m2;
C2 = C1 +m3;
C3 = C2 +m4;
C4 = C3 +m5;
C5 = C4 +m6;




% Obtaining F for Main Tempr Calculations
for i = 2: (C0)
    Fp(i) = 0;
    
end

for i =  (C0)+1:(C1)
    Fp(i)   = 0;
    
end

for i = (C1)+1:(C2)
    Fp(i)   = (-Jhc);
    
end

for i = ((C2)+1):(C3)
    Fp(i)   = (-Jhp);
    
end

for i = ((C3)+1):(C4)
    Fp(i)   = (-Jhbp);
    
end

for i = (C4)+1:(C5)
    Fp(i)   = 0;
    
end

for i = (C5)+1:(m)
    Fp(i)   = 0;
    
end

wp = (b1t + b2t + b3t +t1t +t2t+ t3t+hp)./(m+1);
Fp(1)= -(kt(1)./(wp.^2)).*Th;
Fp(m)= (-kt(end)./(wp.^2)).*Tc;

A = (1./(wp.^2)).*gallery('tridiag', d1,d2,d1);

Tempr = A\Fp';
Tempr(2:end+1)= Tempr;
Tempr(1) = Th;
Tempr(end+1) =Tc;


%% P- Leg Power Output and efficiency

shp = hp./(m4);
shn = hn./(m4);
snm2 = (sn(Tempr(C2+1)) - sn(Tempr(C3)))./(Tempr(C2+1) -Tempr(C3));
snm2 = double(snm2);
spm2 = (sp(Tempr(C2+1)) - sp(Tempr(C3)))./(Tempr(C2+1) -Tempr(C3));
spm2 = double(spm2);

ppm2 =  (pp(Tempr(C2+1)) - pp(Tempr(C3)))./(Tempr(C2+1) -Tempr(C3));  %(pp(Tp(1)) - pp(Tp(Co)))./(Tn(1)-Tn(Co));
ppm2 = double(ppm2);
pnm2 =  (pn(Tempr(C2+1)) - pn(Tempr(C3)))./(Tempr(C2+1) -Tempr(C3));
pnm2 = double(pnm2);

rp2 = (ppm2 .* hp )./(lp.*tp);  % [Ohms] Resistance p-leg
rn2 = (pnm2 .* hn )./(ln.*tn);
S2 = spm2 - snm2;
Rt2 = rp2 + rn2;
RL2 = Rt2 + RC +RBp +RBn; % [Ohms] Load resistance is set to internal resistance
R2 = RL2 +Rt2 +RC+RBp+RBn;
UR =  Rt2 + RC +RBp + RBn;
I2 = (S2.*(Tempr(C2+1) -Tempr(C3)))./(R2);% [A] Current through the unicouple
% I2= NewI;

Voltage = I2*UR;
for j = 1:(m4-1)
    knm(j) = (kn(Tempr(j+C2)) - kn(Tempr(j+1+C2)))./(Tempr(j+C2)-Tempr(j+1+C2));
    
    kpm(j) = (kp(Tempr(j+C2)) - kp(Tempr(j+1+C2)))./(Tempr(j+C2)-Tempr(j+1+C2));
    % kpm(j) = double(kpm);
    snm(j) =  (sn(Tempr(j+C2)) - sn(Tempr(j+1+C2)))./(Tempr(j+C2)-Tempr(j+1+C2));
    % snm = double(snm);
    spm(j) =  (sp(Tempr(j+C2)) - sp(Tempr(j+1+C2)))./(Tempr(j+C2)-Tempr(j+1+C2));
    % spm = double(spm);
    ppm(j) =  (pp(Tempr(j+C2)) - pp(Tempr(j+1+C2)))./(Tempr(j+C2)-Tempr(j+1+C2));
    %ppm = double(ppm);
    pnm(j) =  (pn(Tempr(j+C2)) - pn(Tempr(j+1+C2)))./(Tempr(j+C2)-Tempr(j+1+C2));
    %pnm = double(pnm);
%     shp = hp/(Co);
%     shn = hn/(Co);

    rp1(j) = (ppm(j).* shp )./(lp.*tp);  % [Ohms] Resistance p-leg
    rn1(j) = (pnm(j).* shn )./(ln.*tn); % [Ohms] Resistance n-leg
    cp(j) =  (kpm(j).* lp.*tp)./shp; %  Thermal Conductance of p-leg
    cn(j) = (knm(j) .* ln.*tn)./shn; % Thermal Conductance of n-leg




%      I2 = 4.7;
    %
    Qhp(j) = spm(j).*Tempr(j+C2)*I2 + cp(j).*(Tempr(j+C2)-Tempr(j+C2+1))-(0.5).*(I2.^2).*rp1(j);  % Heat input to hot side
    Qcp(j) = spm(j).*Tempr(j+1+C2).*I2 + cp(j).*(Tempr(j+C2)-Tempr(j+1+C2))+ 0.5.*(I2.^2).*rp1(j); % Heat rejection from cold side
    Pp(j) = (Qhp(j)- Qcp(j)) ;% Power output per leg
   

    Qhn(j) = abs(snm(j)).*Tempr(j+C2).*I2 + cn(j).*(Tempr(j+C2)-Tempr(j+1+C2)) - (0.5).*(I2.^2).*rn1(j) ; % Heat input to hot side
    Qcn(j) = abs(snm(j)).*Tempr(j+1+C2).*I2 + cn(j).*(Tempr(j+C2)-Tempr(j+1+C2)) + 0.5.*(I2.^2).*rn1(j) ;% Heat rejection from cold side
    Pn(j) = Qhn(j)- Qcn(j); % Power output per leg


end
%%


% Joule Heat Top Plate T1
pcm1 =  (pc(Tempr(C1+1)) -pc(Tempr(C2)))./(Tempr(C1+1) -Tempr(C2));
pcm1 = double(pcm1);
rcm1 = (pcm1.*t1w)./(t1l.*t1t);
Jhc1 = I2.^2.*(rcm1);
    
% Joule Heat Bottom Plate B1
pncm2 =  (pc(Tempr(C3+1)) -pc(Tempr(C4)))./(Tempr(C3+1) -Tempr(C4));
pncm2 = double(pncm2);
rncm2 = (pncm2.*b1w)./(b1l.*b1t);
Jhc2 = I2.^2.*(rncm2);



TPp =sum(Pp) -(Jhc1./2)-Jhc2;  % Total Power - P-Leg: Reduce half the Joule Heat from Top Plate from Qh and add the Joule Heat from the bottom to Qc
TPn =sum(Pn) -(Jhc1./2)-Jhc2;% Total Power - N-Leg
%%
Qhin = max(Qhp)+ max(Qhn) -(Jhc1)
P = TPp +TPn %-2.*Jhc2
n =P.*100./Qhin

% figure(1)
% plot(ATemp,'b.')
% hold on
% plot(Tempr,'r.')
% toc

Tempr(m*0.05 +1);
Tempr(m*0.2);
Tempr(C4 +1);
Tempr(C5);

fprintf('Leg Hot Side Temperature is %d C\n', Tempr(m*0.25 +1))
fprintf('Leg Cold Side Temperature is %d C\n', Tempr(m*0.75 ))
fprintf('Current is %d A\n', I2)
fprintf('Voltage is %d V\n', Voltage)

Pv(p) = P;
nv(p) = n;
Qv(p) = Qhin;
Iv(p) = I2;
LH(p) = Tempr(m*0.25 +1);
LC(p) = Tempr(m*0.75 );
Vv(p) = Voltage;
end



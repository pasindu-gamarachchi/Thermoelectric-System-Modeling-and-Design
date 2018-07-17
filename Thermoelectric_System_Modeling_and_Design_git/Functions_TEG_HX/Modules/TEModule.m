% %% Power Generation in a TE Couple, with Power output of 
% two legs seperated
% Only Thermoelectric Legs considered
% Updated on 06/14, material Properties changed to Bi2Te
% Updated on 06/24, material Properties changed to HH, measured by Nick on 05/24
% Properties changed to match that of the paper
% 0718 properties changed to Bi2Te3

  function [QHT, QCT, PT, N,Volt, I, Qm]  = TEModule(Th, Tc, C)

global hn ln tn hp lp tp

% Naming Convention
% property/leg/otherconsideration   no slashes
% hn = 1.69*10^-3;
% % hn = 0.2*(10^-3);
% ln = 1.3*(10^-3);
% tn = 1.3*(10^-3);
% hp = hn;
% lp = 1.3*(10^-3);
% tp = 1.3*(10^-3);
% 
% Th = 300;
% Tc = 160;
% C =1

C_pho = 10*(10^-6); % Contact Resistivity [�?-cm^2] 
%% Temperature dependent properties/ Integrals of functions with respesct to T
%  Resistivity N
syms x
pnf =  int( -6.44051958E-15*x^4 + 3.09525013E-12*x^3 - 4.64686966E-10*x^2 + 5.44900815E-08*x + 3.05911247E-05);
%int( -5.5777E-17*x^5 + 2.8420E-14*x^4 - 5.2132E-12*x^3 + 4.7220E-10*x^2 + 4.9235E-09*x + 3.1567E-05);
%int(-1.89504849E-17*x^4 - 5.64823793E-13*x^3 + 1.48855690E-10*x^2 + 3.25428763E-08*x + 9.14809738E-06); % Original
%int( -6.44051958E-15*x^4 + 3.09525013E-12*x^3 - 4.64686966E-10*x^2 + 5.44900815E-08*x + 3.05911247E-05);
% int(-1.89504849E-17*x^4 - 5.64823793E-13*x^3 + 1.48855690E-10*x^2 + 3.25428763E-08*x + 9.14809738E-06); %int(-1.16561436E-14*x^3 + 5.03594346E-12*x^2 + 6.47260002E-09*x + 6.50670403E-06);%int(-1.16561436E-14*x^3 + 5.03594346E-12*x^2 + 6.47260002E-09*x + 6.50670403E-06);

% Bi2Te3 = int(-1.89504849E-17*x^4 - 5.64823793E-13*x^3 + 1.48855690E-10*x^2 + 3.25428763E-08*x + 9.14809738E-06); %int(-1.16561436E-14*x^3 + 5.03594346E-12*x^2 + 6.47260002E-09*x + 6.50670403E-06);
% Bi2Te3 _film0719 = int( -6.44051958E-15*x^4 + 3.09525013E-12*x^3 - 4.64686966E-10*x^2 + 5.44900815E-08*x + 3.05911247E-05);
% HH = int(3.79548552E-22*x^6 - 6.75630691E-19*x^5 + 4.62762403E-16*x^4 - 1.59187729E-13*x^3 + 2.48272530E-11*x^2 + 4.34512734E-09*x + 6.47172345E-06);
%% Seebeck N
syms x
snf =  int( 4.02004040E-14*x^4 - 2.01825724E-11*x^3 + 4.01229697E-09*x^2 - 4.47528114E-07*x - 1.16757881E-04);
%int ( 4.0200E-08*x^4 - 2.0183E-05*x^3 + 4.0123E-03*x^2 - 4.4753E-01*x - 1.1676E+02);
%int(3.28019971E-14*x^4 - 1.66155332E-11*x^3 + 4.62835101E-09*x^2 - 6.11524875E-07*x - 1.76396851E-04);  % Original
%int(5.48885661E-14*x^3 + 2.22674400E-10*x^2 - 2.66412731E-07*x - 1.33513735E-04);int(-3.07307824E-20*x^6 + 5.17273632E-17*x^5 - 3.14308248E-14*x^4 + 8.23370332E-12*x^3 - 6.95427569E-10*x^2 - 2.15931357E-07*x - 1.30728621E-04);
%Bi2Te3 =int(3.28019971E-14*x^4 - 1.66155332E-11*x^3 + 4.62835101E-09*x^2 - 6.11524875E-07*x - 1.76396851E-04); %int(5.48885661E-14*x^3 + 2.22674400E-10*x^2 - 2.66412731E-07*x - 1.33513735E-04);
% Bi2Te3 _film0719 = int( 4.02004040E-14*x^4 - 2.01825724E-11*x^3 + 4.01229697E-09*x^2 - 4.47528114E-07*x - 1.16757881E-04)

% HH = int(-3.07307824E-20*x^6 + 5.17273632E-17*x^5 - 3.14308248E-14*x^4 + 8.23370332E-12*x^3 - 6.95427569E-10*x^2 - 2.15931357E-07*x - 1.30728621E-04);
%% Conductivity N
syms x
knf =  int( 6.77463453E-11*x^4 - 3.55085699E-08*x^3 + 1.16901257E-05*x^2 - 1.39832195E-03*x + 5.90438185E-01);
%int ( 6.7746E-11*x^4 - 3.5509E-08*x^3 + 1.1690E-05*x^2 - 1.3983E-03*x + 5.9044E-01); 
%int(-1.25842544E-10*x^4 + 8.58521996E-08*x^3 - 5.60764755E-06*x^2 - 8.78714664E-04*x + 1.28916143E+00); % Original
 % int(1.18722952E-08*x^3 - 4.65123739E-06*x^2 - 7.31176479E-04*x + 3.53715861E+00);int( -3.42296835E-16*x^6 + 5.14867564E-13*x^5 - 1.96473071E-10*x^4 - 4.33833033E-08*x^3 + 4.91278625E-05*x^2 - 1.38542448E-02*x + 5.56989213E+00);
%Bi2Te3 =int(-1.25842544E-10*x^4 + 8.58521996E-08*x^3 - 5.60764755E-06*x^2 - 8.78714664E-04*x + 1.28916143E+00); % int(1.18722952E-08*x^3 - 4.65123739E-06*x^2 - 7.31176479E-04*x + 3.53715861E+00);
% Bi2Te3 _film0719 = int( 6.77463453E-11*x^4 - 3.55085699E-08*x^3 + 1.16901257E-05*x^2 - 1.39832195E-03*x + 5.90438185E-01);
% HH = int( -3.42296835E-16*x^6 + 5.14867564E-13*x^5 - 1.96473071E-10*x^4 - 4.33833033E-08*x^3 + 4.91278625E-05*x^2 - 1.38542448E-02*x + 5.56989213E+00);
%% Resisitivity P
syms x
ppf = int ( -3.61477E-15*x^4 + 3.09506E-12*x^3 - 7.43539E-10*x^2 + 1.16363E-07*x + 1.72759E-05);
%int( -1.28936729E-15*x^4 + 1.18635837E-14*x^3 + 1.38861195E-10*x^2 + 4.81304310E-08*x + 8.63224270E-06); % Original
% int(-1.26642460E-15*x^3 - 6.37965919E-12*x^2 + 1.55655066E-08*x + 1.03160255E-05);
% 
% int(-1.26642460E-15*x^3 - 6.37965919E-12*x^2 + 1.55655066E-08*x + 1.03160255E-05);
%int(1.26762914E-22*x^6 - 2.35719632E-19*x^5 + 1.60899624E-16*x^4 - 5.53161599E-14*x^3 + 1.53924136E-11*x^2 + 4.75622122E-09*x + 1.87165436E-06); 
% int(-1.26642460E-15*x^3 - 6.37965919E-12*x^2 + 1.55655066E-08*x + 1.03160255E-05);
%Bi2Te3 =int( -1.28936729E-15*x^4 + 1.18635837E-14*x^3 + 1.38861195E-10*x^2 + 4.81304310E-08*x + 8.63224270E-06); % int(-1.26642460E-15*x^3 - 6.37965919E-12*x^2 + 1.55655066E-08*x + 1.03160255E-05);
% Bi2Te3 _film0719 = int( -3.61477286E-15*x^4 + 3.09505543E-12*x^3 - 7.43538765E-10*x^2 + 1.16362916E-07*x + 1.72758675E-05); 
% HH = int(1.26762914E-22*x^6 - 2.35719632E-19*x^5 + 1.60899624E-16*x^4 - 5.53161599E-14*x^3 + 1.53924136E-11*x^2 + 4.75622122E-09*x + 1.87165436E-06);
%% Seebeck P
 syms x
 spf = int( -1.58933333E-14*x^4 + 1.23395556E-11*x^3 - 4.23451111E-09*x^2 + 7.13504603E-07*x + 1.61110500E-04);
 %int ( -1.58933E-14*x^4 + 1.23396E-11*x^3 - 4.23451E-09*x^2 + 7.13505E-07*x + 1.61110E-04);
 %int( -2.16011384E-14*x^4 + 1.21036642E-11*x^3 - 4.24902245E-09*x^2 + 7.16773656E-07*x + 1.75862621E-04); % Original
% 
 %int( -2.16011384E-14*x^4 + 1.21036642E-11*x^3 - 4.24902245E-09*x^2 + 7.16773656E-07*x + 1.75862621E-04); %int( 1.50803298E-13*x^3 - 3.37285775E-10*x^2 + 3.01002672E-07*x + 1.23997961E-04);int(-1.36579995E-20*x^6 + 2.24957017E-17*x^5 - 1.38487487E-14*x^4 + 3.96125625E-12*x^3 - 5.38600018E-10*x^2 + 2.14193134E-07*x + 6.98593256E-05);
 %Bi2Te3 =int( -2.16011384E-14*x^4 + 1.21036642E-11*x^3 - 4.24902245E-09*x^2 + 7.16773656E-07*x + 1.75862621E-04); %int( 1.50803298E-13*x^3 - 3.37285775E-10*x^2 + 3.01002672E-07*x + 1.23997961E-04);
% Bi2Te3 _film0719 = int( -1.58933333E-14*x^4 + 1.23395556E-11*x^3 - 4.23451111E-09*x^2 + 7.13504603E-07*x + 1.61110500E-04);
 % HH = int(-1.36579995E-20*x^6 + 2.24957017E-17*x^5 - 1.38487487E-14*x^4 + 3.96125625E-12*x^3 - 5.38600018E-10*x^2 + 2.14193134E-07*x + 6.98593256E-05);
 %% Conductivity P
 
 syms x
 kpf =  int( -7.34772339E-14*x^6 + 7.14377389E-11*x^5 - 2.71158280E-08*x^4 + 5.05485382E-06*x^3 - 4.64386257E-04*x^2 + 1.82440425E-02*x + 5.85287146E-01);
 %int ( -2.93820E-10*x^4 + 5.48637E-08*x^3 + 2.33659E-05*x^2 - 5.27538E-03*x + 1.02088E+00);
 %int( -1.92004259E-10*x^4 + 1.19050111E-07*x^3 - 1.45459874E-05*x^2 - 4.95909112E-04*x + 1.10059563E+00); %Original
% 

 %Bi2Te3 =int( -1.92004259E-10*x^4 + 1.19050111E-07*x^3 - 1.45459874E-05*x^2 - 4.95909112E-04*x + 1.10059563E+00);%int( -1.79754155E-16*x^6 + 3.96173694E-13*x^5 - 3.42301654E-10*x^4 + 1.48422824E-07*x^3 - 3.25292137E-05*x^2 + 2.66399769E-03*x + 3.06261874E+00);
 
 % Bi2Te3 _film0719 = int( -7.34772339E-14*x^6 + 7.14377389E-11*x^5 - 2.71158280E-08*x^4 + 5.05485382E-06*x^3 - 4.64386257E-04*x^2 + 1.82440425E-02*x + 5.85287146E-01);
 %HH = int( -8.25275365E-16*x^6 + 1.50160725E-12*x^5 - 1.02862881E-09*x^4 + 3.18113901E-07*x^3 - 3.14899647E-05*x^2 - 9.18620390E-03*x + 7.73547497E+00);
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
knm =  double(knm); %  0.77 ;
kpm = (kp(Th) - kp (Tc))./(Th-Tc);
kpm =   double(kpm); % 0.77 ;
snm =  (sn(Th) - sn(Tc))./(Th-Tc);
snm =  double(snm); % -133.5*(10^-6);
spm =  (sp(Th) - sp(Tc))./(Th-Tc);
spm =  double(spm); % 133.5*(10^-6) ;%
ppm =  (pp(Th) - pp(Tc))./(Th-Tc);
ppm =   double(ppm); % 20*(10^-6);
pnm =  (pn(Th) - pn(Tc))./(Th-Tc);
pnm =   double(pnm); % 20*(10^-6);

%%

% wp = hp./(m+1); % mesh height p-leg
% e =ones(m+2,1);
% A = (1/wp.^2).*spdiags( [e -2.*e e], -1:1, m, m);
% Fp = (-1./kpm(i)).*ones(m,1);
% Fp(1)= (-1./kpm(i)) - Th(i)./(wp.^2);
% Fp(m)= (-1./kpm(i)) - Tc./(wp.^2);

% Tp =A\Fp;
% 
% wn = hn./(m+1);
% Fn = (-1./knm(i)).*ones(m,1);
% Fn(1)= (-1./knm(i)) - Th(i)./(wn.^2);
% Fn(m)= (-1./knm(i)) - Tc./(wn.^2);
% 
% z = linspace(hn -wn,0+wn,m);
% 
% Tn =A\Fn;

%% Calculations

rp = (ppm .* hp )./(lp.*tp);  % [Ohms] Resistance p-leg
rn = (pnm .* hn )./(ln.*tn); % [Ohms] Resistance n-leg
cp =  (kpm .* lp.*tp)./hp; %  Thermal Conductance of p-leg
cn = (knm .* ln.*tn)./hn; % Thermal Conductance of n-leg

CR_n = C_pho/((tn*100)*(ln*100));
CR_p = C_pho/((tp*100)*(lp*100));

CR = 2*(CR_n + CR_p);

% CR = (C_pho*(ln*100)*(tn*100) + C_pho*(lp*100)*(tp*100))*2;
Rt = rp +rn + CR; % [Ohms] Resistance Total 
RL = Rt; % [Ohms] Load resistance is set to internal resistance
R = RL +Rt; % [Ohms] 
S =  spm -snm ; % Seebeck of couple
I = (S.*(Th -Tc))./(R);% [A] Current through the unicouple
V = I*Rt;
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
N = P/(Qhp +Qhn)*100;
%  end
Qh = Qhp + Qhn;
Qc = Qcp + Qcn;
Pn;
Pp;
P;

QHT = Qh*C;
QCT = Qc*C;
PT = P*C;
Volt = V*C;
Qm = (QHT +QCT)/2;

 end
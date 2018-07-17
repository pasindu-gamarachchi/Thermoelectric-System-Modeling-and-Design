% Cold Side Heat Exchanger
% Function Calculates Heat Transferred from 
% cold side fins to cold fluid 


% Working fluid changed to air

% Updated on 04/28
% Replaced Nickel with Aluminum Pure
  function [Qc, h, TEMC, Ra] = HXFreeConv(Tcb, Tamb)

global NC tf_C  cvl Tcg Tcbg hxlC hxhC hxwC  Cb

g = 9.81;

syms x

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

kAlpf = int( 7.00000000E-09*x^4 - 6.18933333E-06*x^3 + 1.25071800E-03*x^2 + 1.47795093E-02*x + 2.28807284E+02);
kAlp = symfun(kAlpf,x);

kCupf = int(5.25E-15*x^6 - 1.34E-11*x^5 + 1.26E-08*x^4 - 5.19E-06*x^3 + 8.71E-04*x^2 - 9.61E-02*x + 3.99E+02);
kCup = symfun(kCupf,x);

Tb = Tamb + Tcbg;

Pf = 2*(cvl);
Af = cvl*tf_C;

kAlpm = (kAlp(Tcb) - kAlp(Tcb-0.5))./(0.5);
kAlpm = double(kAlpm);

kCupm = (kCup(Tcb) - kCup(Tcb-0.5))./(0.5);
kCupm = double(kCupm);
% Set operating Cold Side Metal
k =  kCupm;


%% Start of Cold convection function
Tcold = Tamb;
Tco = Tcold + Tcg;
%%  Temperature dependent properties

 syms x  
 % Density of water/ Hot side fluid

 daf = int(2.876602E-12*x^4 - 7.350893E-09*x^3 + 7.284062E-06*x^2 - 3.760334E-03*x + 1.251051E+00);
% Specific Heat Capacity of water
%  cwf = int( 5.951927E-06*x^4 - 1.221571E-03*x^3 + 9.793076E-02*x^2 - 3.302581E+00*x + 4.217176E+03);
 caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03);
% Dynamic viscocity of Water
%  mwf = int(4.916802E-11*x^4 - 1.125020E-08*x^3 + 1.050048E-06*x^2 - 5.463958E-05*x + 1.756519E-03);
 maf = int(-1.17230617E-17*x^4 +3.14251436E-14*x^3 -3.87294440E-11*x^2 +4.96283182E-08*x + 1.71301016E-05); 
% Thermal Cond of Water
%  kwf = int( 1.934482E-10*x^4 - 4.145917E-08*x^3 - 4.292555E-06*x^2 + 1.765088E-03*x + 5.695880E-01);
kaf =  int( 1.146049E-14*x^4 + 5.825050E-12*x^3 - 4.412201E-08*x^2 + 8.330657E-05*x + 2.396400E-02);
% Prandtl Number of Water
%  pwf = int( 4.299125E-07*x^4 - 9.796954E-05*x^3 + 9.010003E-03*x^2 - 4.486901E-01*x + 1.305137E+01);
 paf = int( 1.946038E-13*x^4 - 8.530084E-10*x^3 + 9.958116E-07*x^2 - 3.351574E-04*x + 7.176453E-01);
% Thermal Diffusivity of Air
faf = int( 1.093417E-16*x^4 - 1.924826E-13*x^3 + 1.754598E-10*x^2 + 1.375914E-07*x + 1.847226E-05);

% Kinematic Viscocity of Air
vaf = int( 9.471589E-18*x^4 - 3.915231E-14*x^3 + 1.054041E-10*x^2 + 8.961175E-08*x + 1.340998E-05);


da =  symfun(daf,x);
ca =  symfun(caf,x);
ma =  symfun(maf,x);
ka =  symfun(kaf,x);
pa =  symfun(paf,x);
fa =  symfun(faf,x); 
va =  symfun(vaf,x);




Tfilm = (Tamb + Tcb)/2;

% Air
dam = (da(Tfilm) - da (Tfilm-0.5))./(0.5);
dam = double(dam);

cam = (ca(Tfilm) - ca(Tfilm-0.5))./(0.5);
cam = double(cam);

mam = (ma(Tfilm) - ma(Tfilm-0.5))./(0.5);
mam = double(mam);

kam = (ka(Tfilm) - ka(Tfilm-0.5))./(0.5);
kam = double(kam);

pam = (pa(Tfilm) - pa(Tfilm-0.5))./(0.5);
pam = double(pam);

fam = (fa(Tfilm) - fa(Tfilm-0.5))./(0.5);
fam = double(fam);

vam = (va(Tfilm) - va(Tfilm-0.5))./(0.5);
vam = double(vam);

TfilmK = Tfilm +273.15;
Be = 1/TfilmK;

% Intermediate Calculations
sC = (hxwC - (NC.*tf_C))./(NC-1);  % Fin Spacing

al = sC./hxlC;


Ra = g*Be*(Tb -Tamb)*(sC^3)/(fam*vam);

Nu =  (576/((Ra*al)^2) + 2.87/((Ra*al)^0.5))^(-0.5);%(1/24)*Ra*(al)*(1 - exp(-35/(Ra*(al))))^(0.75); % 

h = (Nu*kam)/sC;
%  h =30;
 M = (sqrt(h*Pf*k*Af))*(Tcb -Tamb);
m = sqrt((h*Pf)/(k*Af));

Qc = M*tanh(m*hxhC)*(NC);

% Tcb_c = Qc/((sqrt(h*Pf*k*Af))*tanh(m*hxhC)*(NC-1)) +Tamb



Rb = Cb/(k*cvl*hxwC);

TEMC = Rb*Qc + Tcb;

%  Tcb =30;
%  Tamb = 20;
% Tcb =70;
% Tcf =20;

% All cold side Heat Exchanger Parameters in this file
% Tcg  = 8;
% 
% 
% Tcbg= 50;
% hxlC = 0.16; % [m]

% hxwC = 0.04; % [m]

% PackFraction = (tf_C*NC/hxwC)

%%
% sC = (hxw - (N_C.*tf_C))./(N_C+1); % Fin Spacing
% 
% cwf = int(5.951927E-06*x^4 - 1.221571E-03*x^3 + 9.793076E-02*x^2 - 3.302581E+00*x + 4.217176E+03); % Specific Heat of Water
% cw = symfun(cwf, x);
% cwm = (cw(Tcg) - cw(Tcf))./(Tcg-Tcf);
% cwm = double(cwm);
% 
% cp =cwm;
% knim = (kni(Tb) - kni(Tb-10))./(10);
% knim = double(knim);
% dw =  symfun(dwf,x);
% cw =  symfun(cwf,x);
% mw =  symfun(mwf,x);
% kw =  symfun(kwf,x);
% pw =  symfun(pwf,x);
% dwm = (dw(Tco) - dw (Tcold))./(Tcg);
% dwm = double(dwm);
% 
% cwm = (cw(Tco) - cw(Tcold))./(Tcg);
% cwm = double(cwm);
% 
% mwm = (mw(Tco) - mw(Tcold))./(Tcg);
% mwm = double(mwm);
% 
% kwm = (kw(Tco) - kw(Tcold))./(Tcg);
% kwm = double(kwm);
% 
% pwm = (pw(Tco) - pw(Tcold))./(Tcg);
% pwm = double(pwm);

% dhC = (4*sC.*(hxhC))./(sC + 2*hxhC );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
% v = mfC./(dwm.*sC.*(hxhC)*(NC-1));
% Re = v.*dwm.*dhC/(mwm);
% Pr = pwm;


% % Colburn and pressure factors
% al = sC./(hxhC);
% if Re<=1000
% j_f = 0.483*(cvl/dhC)^(-0.162)*al^(-0.184)*Re^(-0.536);
% f_f = 7.661*(cvl/dhC)^(-0.384)*al^(-0.092)*Re^(-0.712);
% 
% elseif Re>=1000
% % j_f = 0.242*(1/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% % f_f = 1.136*(l/Dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% % j_f = 0.242*(cvl/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% 
% j_f =0.5*(cvl/dhC)^(-0.322)*(tf_C/dhC)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
% f_f = 1.136*(cvl/dhC)^(-0.3)*(tf_C/dhC)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
% end
% 
% 
% 
%  h = j_f.*Re.*(Pr.^(1./3)).*kwm./(dhC);
% Tcb
% Tcf
% cwm
% mfC
% Tout = Tamb + Qc/(mfC*cwm);
% h


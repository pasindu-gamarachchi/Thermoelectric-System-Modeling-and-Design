% Natural Convection Cold Side Heat Exchanger
% Function Calculates Heat Transferred from Heat Sink Fins
% Calculations based on source : http://www.electronics-cooling.com/2002/02/estimating-natural-convection-heat-transfer-for-arrays-of-vertical-parallel-flat-plates/


  function [Qc h TEMC Ra] = HXFreeConvV3(Tcb, Tamb)

global NC tf_C  cvl  hxlC hxhC hxwC  Cb

g = 9.81;

syms x

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

kAlpf = int( 7.00000000E-09*x^4 - 6.18933333E-06*x^3 + 1.25071800E-03*x^2 + 1.47795093E-02*x + 2.28807284E+02);
kAlp = symfun(kAlpf,x);

kCupf = int(5.25E-15*x^6 - 1.34E-11*x^5 + 1.26E-08*x^4 - 5.19E-06*x^3 + 8.71E-04*x^2 - 9.61E-02*x + 3.99E+02);
kCup = symfun(kCupf,x);


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





kAlpm = (kAlp(Tcb) - kAlp(Tcb-0.5))./(0.5);
kAlpm = double(kAlpm);

kCupm = (kCup(Tcb) - kCup(Tcb-0.5))./(0.5);
kCupm = double(kCupm);
% Set operating Cold Side Metal
k =  kCupm;


%% Start of Cold convection function



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

Pf = 2*(cvl);
Af = cvl*tf_C;


% Intermediate Calculations
sC = (hxwC - (NC.*tf_C))./(NC-1);  % Fin Spacing

if sC <= 0
    fprintf('Spacing less than 1\n')
    return
end

al = sC./hxlC;


Ra = ((dam)^2)*g*Be*cam*(sC^4)*(Tcb-Tamb)/(mam*kam*hxlC);

Nu = (Ra/24)*(1 -exp(-35/Ra))^(0.75);

h = (Nu*kam)/sC;
%  h =30;

Qc = h*(NC*2*hxhC*hxlC)*(Tcb-Tamb);


% Tcb_c = Qc/((sqrt(h*Pf*k*Af))*tanh(m*hxhC)*(NC-1)) +Tamb



Rb = Cb/(k*cvl*hxwC);

TEMC = Rb*Qc + Tcb;


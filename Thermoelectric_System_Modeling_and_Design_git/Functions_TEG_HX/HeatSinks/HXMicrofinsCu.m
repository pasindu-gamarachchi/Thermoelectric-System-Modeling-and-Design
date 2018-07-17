% Cold Side Heat Exchanger
% Function Calculates Heat Transferred from 
% cold side fins to cold fluid 

% Microfins, with enhanced convection used
% 06/01 spacing equation updated
% Convection coefficients from source used : [1] Guan, N., Liu, Z., Zhang, C., and Jiang, G., 2013, "Natural convection heat transfer on surfaces of copper micro-wires", Heat and Mass Transfer, 50(2), pp. 275-284.



  function [Qc, h ,TEMC, Gr ] = HXMicrofinsCu(Tcb, Tamb)

global  cvl Tcg hxwC  Cb NCx NCy fd fh

% Tcb = 26;
% Tamb = 22;
% 
% fd = 57*(10^-6);
% cvl =0.04;
% Tcg = 0.01;
% hxwC = 0.04;
% Cb  = 5*(10^-3);
% NCx = 30;
% NCy = 30;
% fh = 5*(10^-3);

g = 9.81;
lam = 68*10^-9;


% Fin Material
syms x

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

kAlpf = int( 7.00000000E-09*x^4 - 6.18933333E-06*x^3 + 1.25071800E-03*x^2 + 1.47795093E-02*x + 2.28807284E+02);
kAlp = symfun(kAlpf,x);

kCupf = int(5.25E-15*x^6 - 1.34E-11*x^5 + 1.26E-08*x^4 - 5.19E-06*x^3 + 8.71E-04*x^2 - 9.61E-02*x + 3.99E+02);
kCup = symfun(kCupf,x);

% Tb = Tamb + Tcbg;

Pf = pi*fd;
Af = (pi/4)*(fd^2);

kAlpm = (kAlp(Tcb) - kAlp(Tcb-0.5))./(0.5);
kAlpm = double(kAlpm);

kCupm = (kCup(Tcb) - kCup(Tcb-0.5))./(0.5);
kCupm = double(kCupm);

% Set operating Cold Side Metal
k = kCupm;


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
% Thermal Cond of Air
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
dam = (da(Tfilm) - da (Tfilm-0.1))./(0.1);
dam = double(dam);

cam = (ca(Tfilm) - ca(Tfilm-0.1))./(0.1);
cam = double(cam);

mam = (ma(Tfilm) - ma(Tfilm-0.1))./(0.1);
mam = double(mam);

kam = (ka(Tfilm) - ka(Tfilm-0.1))./(0.1);
kam = double(kam);

pam = (pa(Tfilm) - pa(Tfilm-0.1))./(0.1);
pam = double(pam);

fam = (fa(Tfilm) - fa(Tfilm-0.1))./(0.1);
fam = double(fam);

vam = (va(Tfilm) - va(Tfilm-0.1))./(0.1);
vam = double(vam);


TfilmK = Tfilm +273.15;
Be = 1/TfilmK;

% Intermediate Calculations
sC = (hxwC -(fd*NCx))/(NCx -1);
% al = sC./hxlC;
% Ra = g*Be*(Tb -Tamb)*(sC^3)/(fam*vam);

if sC < 0.0001
    fprintf('Spacing Condition not satisfied\n')
    return
else
end


Gr = g*Be*(Tcb- Tamb)*(fd^3)/(vam^2);

if Gr < 0.0001 | Gr > 2.5
    fprintf('Gr Condition not satisfied\n')
    Gr
%     return
else
end

Nu = 1.03*(Gr*pam)^(0.035);

h = Nu*(kam)/(fd);


M = (sqrt(h*Pf*k*Af))*(Tcb -Tamb);
m = sqrt((h*Pf)/(k*Af));

Qc = M*tanh(m*fh)*(NCx*NCy);

% Tcb_c = Qc/((sqrt(h*Pf*k*Af))*tanh(m*hxhC)*(NC-1)) +Tamb



Rb = Cb/(k*cvl*hxwC);

TEMC = Rb*Qc + Tcb;




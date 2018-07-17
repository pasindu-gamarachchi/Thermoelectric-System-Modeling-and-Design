% Natural Convection Cold Side Heat Exchanger
% Function Calculates Heat Transferred from Heat Sink Fins
% Calculations based on source : [1] Bar-Cohen, A., Iyengar, M., and Kraus, A., 2003, "Design of Optimum Plate-Fin Natural Convective Heat Sinks", Journal of Electronic Packaging, 125(2), p. 208.
% Function Calculates Total Heat Tranfer, considering both the base and the
% Heat transfer from Fins, also calculates the Tip Temperature
% Jun 23, Radiation effects added, Average Fin Temp used for Radiation Hot
% Temperature

  function [QT, Q_b, h_f, h_b, TEMC, TTip, et_f, Qr, sC] = HXFreeConvV4(Tcb, Tamb)


global NC tf_C  cvl  hxlC  hxwC  Cb hxhC

g = 9.81;
Boltz = 5.67*10^-8;


syms x

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

kAlpf = int( 7.00000000E-09*x^4 - 6.18933333E-06*x^3 + 1.25071800E-03*x^2 + 1.47795093E-02*x + 2.28807284E+02);
kAlp = symfun(kAlpf,x);

kCupf = int(5.25E-15*x^6 - 1.34E-11*x^5 + 1.26E-08*x^4 - 5.19E-06*x^3 + 8.71E-04*x^2 - 9.61E-02*x + 3.99E+02);
kCup = symfun(kCupf,x);

ep_Cu =0.28;

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



% Intermediate Calculations
sC = (hxwC - (NC.*tf_C))./(NC-1); % Fin Spacing

if sC <= 0
    fprintf('Spacing less than 0\n')
    return
end
al = sC./hxlC;


Ra = g*Be*(Tcb -Tamb)*(sC^3)/(fam*vam);


Ra_b = (Tcb-Tamb)*g*Be*(pam)*(hxlC^3)/(vam^2); % Raleigh Number for Base Area

Nu_b = 0.59*(Ra_b^(0.25)); % Nusselt Number for Base Area

h_b = (Nu_b*kam)/(hxlC); % Convection Coeff from Base Area 

Q_b = (sC*hxlC*(NC-1))*(h_b)*(Tcb-Tamb); % Heat transfer from Base Area

% El = g*Be*(Tcb-Tamb)*(pam)*(sC^4)/(hxlC*vam^2) % Elenbaas Number

et_fin = 0.01; % Fin Efficiency

Nu_f =  (576/((et_fin*Ra*al)^2) + 2.87/((et_fin*Ra*al)^0.5))^(-0.5);

h_f = Nu_f*kam/(sC);

m = sqrt(2*h_f/(k*tf_C));

et_f = tanh(m*hxhC)/(m*hxhC);
%%
while abs(et_fin -et_f) > 0.001
    et_fin = et_fin +0.001;
    Nu_f =  (576/((et_fin*Ra*al)^2) + 2.87/((et_fin*Ra*al)^0.5))^(-0.5);%(1/24)*Ra*(al)*(1 - exp(-35/(Ra*(al))))^(0.75); % 

    h_f = Nu_f*kam/(sC);

    m = sqrt(2*h_f/(k*tf_C));

    et_f = tanh(m*hxhC)/(m*hxhC);
    
    if et_fin >1
        fprintf('Error!\n')
        return
    end

end

Q_f = hxlC*k*tf_C*(Tcb - Tamb)*m*(tanh(m*hxhC));  % Similar to Qf = Mtanh(mL)
% Nu =  (576/((Ra)^2) + 2.873/((Ra)^0.5))^(-0.5)%(1/24)*Ra*(al)*(1 - exp(-35/(Ra*(al))))^(0.75); % 

TTip = (Tcb-Tamb)/(cosh(m*hxhC)) + Tamb;

% Radiation 
A = (hxhC*hxlC*2) + (hxhC*tf_C*NC*2) + (hxlC*tf_C*NC) +(Cb*hxwC*2)+(Cb*hxlC*2);

Tamb_K = Tamb + 273.15;
T_Rad = mean([ Tcb, TTip]);
T_RadK = T_Rad + 273.15;

Qr = A*(Boltz)*(ep_Cu)*(T_RadK^4 - Tamb_K^4);



QT = Q_f*NC + Q_b +Qr;


Rb = Cb/(k*cvl*hxwC);

TEMC = Rb*QT + Tcb;



% Microfins convection coefficient

  function h = convectionmicro(Tcb, Tamb)

global  fd


g = 9.81;
lam = 68*10^-9;


% Fin Material
syms x

knif = ( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

kAlpf = ( 7.00000000E-09*x^4 - 6.18933333E-06*x^3 + 1.25071800E-03*x^2 + 1.47795093E-02*x + 2.28807284E+02);
kAlp = symfun(kAlpf,x);

kCupf = (5.25E-15*x^6 - 1.34E-11*x^5 + 1.26E-08*x^4 - 5.19E-06*x^3 + 8.71E-04*x^2 - 9.61E-02*x + 3.99E+02);
kCup = symfun(kCupf,x);

kAlpm = kAlp(Tcb);
kAlpm = double(kAlpm);

kCupm = kCup(Tcb) ;
kCupm = double(kCupm);


k = kCupm;



Tcold = Tamb;

%%  Temperature dependent properties

 syms x  
 daf = (2.876602E-12*x^4 - 7.350893E-09*x^3 + 7.284062E-06*x^2 - 3.760334E-03*x + 1.251051E+00);

 caf = ( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03);

 maf = (-1.17230617E-17*x^4 +3.14251436E-14*x^3 -3.87294440E-11*x^2 +4.96283182E-08*x + 1.71301016E-05); 

kaf =  ( 1.146049E-14*x^4 + 5.825050E-12*x^3 - 4.412201E-08*x^2 + 8.330657E-05*x + 2.396400E-02);

 paf = ( 1.946038E-13*x^4 - 8.530084E-10*x^3 + 9.958116E-07*x^2 - 3.351574E-04*x + 7.176453E-01);

faf = ( 1.093417E-16*x^4 - 1.924826E-13*x^3 + 1.754598E-10*x^2 + 1.375914E-07*x + 1.847226E-05);

% Kinematic Viscocity of Air
vaf = ( 9.471589E-18*x^4 - 3.915231E-14*x^3 + 1.054041E-10*x^2 + 8.961175E-08*x + 1.340998E-05);

da =  symfun(daf,x);
ca =  symfun(caf,x);
ma =  symfun(maf,x);
ka =  symfun(kaf,x);
pa =  symfun(paf,x);
fa =  symfun(faf,x); 
va =  symfun(vaf,x);


Tfilm = (Tamb + Tcb)/2;

% Air
dam = da(Tfilm) ;
dam = double(dam);

cam = ca(Tfilm);
cam = double(cam);

mam = ma(Tfilm) ;
mam = double(mam);

kam = ka(Tfilm) ;
kam = double(kam);

pam = pa(Tfilm) ;
pam = double(pam);

fam = fa(Tfilm) ;
fam = double(fam);

vam = va(Tfilm);
vam = double(vam);


TfilmK = Tfilm +273.15;
Be = 1/TfilmK;


% Intermediate Calculations
Kn = lam/fd;
C = fam/(fd^2);
C1 = (log(C))^2;
C2 = log(C);


h = (kam/(1+Kn))*(1/fd)*((1/16)*C1 -0.292*C2 +0.958)^(-0.5);

  end



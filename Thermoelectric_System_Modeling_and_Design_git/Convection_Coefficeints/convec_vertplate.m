% Convection Coefficient for heat transfer from upward facing hot plate
% Date Created : 02/22/2017

 function h_b = convec_vertplate(Tcb ,Tamb, Char_L)

g = 9.81;

 syms x  

% Thermal Cond of Air
%  kwf = int( 1.934482E-10*x^4 - 4.145917E-08*x^3 - 4.292555E-06*x^2 + 1.765088E-03*x + 5.695880E-01);
kaf =  ( 1.146049E-14*x^4 + 5.825050E-12*x^3 - 4.412201E-08*x^2 + 8.330657E-05*x + 2.396400E-02);
% Prandtl Number of Water
%  pwf = int( 4.299125E-07*x^4 - 9.796954E-05*x^3 + 9.010003E-03*x^2 - 4.486901E-01*x + 1.305137E+01);
 paf = ( 1.946038E-13*x^4 - 8.530084E-10*x^3 + 9.958116E-07*x^2 - 3.351574E-04*x + 7.176453E-01);
% Kinematic Viscocity of Air
vaf = ( 9.471589E-18*x^4 - 3.915231E-14*x^3 + 1.054041E-10*x^2 + 8.961175E-08*x + 1.340998E-05);



ka =  symfun(kaf,x);
pa =  symfun(paf,x);
va =  symfun(vaf,x);

Tfilm = (Tamb + Tcb)/2;

% Air

kam = ka(Tfilm);
kam = double(kam);

pam = pa(Tfilm) ;
pam = double(pam);

vam = va(Tfilm) ;
vam = double(vam);

TfilmK = Tfilm +273.15;
Be = 1/TfilmK;

Ra_b = (Tcb-Tamb)*g*Be*(pam)*(Char_L^3)/(vam^2); % Raleigh Number for Base Area

Nu_b = 0.59*(Ra_b^(0.25)); % Nusselt Number for Base Area

h_b = (Nu_b*kam)/(Char_L); % Convection Coeff from Base Area 

 end
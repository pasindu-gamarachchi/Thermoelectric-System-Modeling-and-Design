% Function calculates the convection coefficient along with Velocity and
% Reynolds Number for the hot side of the heat exchanger which uses air as
% the working fluid 
function [h, v, Re] = coldconvcoeff(Tcold, mfC)

%  global hxhC hxwC NC tf_C   cvl Tcg
% Tcold =5;
% mfC = 0.0005;
 


% 
%  hxl =0.04;
  hxhC =0.2;
 hxwC= 0.04;
 NC=5;
 tf_C=0.1*10^-3;
  cvl =3*10^-3;
  Tcg =6;


Tco = Tcold + Tcg;
%%  Temperature dependent properties

 syms x  
 % Density of water/ Hot side fluid
dwf = int(-1.244106E-07*x^4 + 4.224596E-05*x^3 - 7.670379E-03*x^2 + 5.169102E-02*x + 1.000067E+03);

% Specific Heat Capacity of water
cwf = int( 5.951927E-06*x^4 - 1.221571E-03*x^3 + 9.793076E-02*x^2 - 3.302581E+00*x + 4.217176E+03);
% Dynamic viscocity of Water
mwf = int(4.916802E-11*x^4 - 1.125020E-08*x^3 + 1.050048E-06*x^2 - 5.463958E-05*x + 1.756519E-03);
% Thermal Cond of Water
kwf = int( 1.934482E-10*x^4 - 4.145917E-08*x^3 - 4.292555E-06*x^2 + 1.765088E-03*x + 5.695880E-01);

% Prandtl Number of Water
pwf = int( 4.299125E-07*x^4 - 9.796954E-05*x^3 + 9.010003E-03*x^2 - 4.486901E-01*x + 1.305137E+01);


dw =  symfun(dwf,x);
cw =  symfun(cwf,x);
mw =  symfun(mwf,x);
kw =  symfun(kwf,x);
pw =  symfun(pwf,x);




% Integral averages

dwm = (dw(Tco) - dw (Tcold))./(Tcg);
dwm = double(dwm);

cwm = (cw(Tco) - cw(Tcold))./(Tcg);
cwm = double(cwm);

mwm = (mw(Tco) - mw(Tcold))./(Tcg);
mwm = double(mwm);

kwm = (kw(Tco) - kw(Tcold))./(Tcg);
kwm = double(kwm);

pwm = (pw(Tco) - pw(Tcold))./(Tcg);
pwm = double(pwm);







% Intermediate Calculations
sC = (hxwC - (NC.*tf_C))./(NC+1);  % Fin Spacing
dhC = (4*sC.*(hxhC))./(sC + 2*hxhC );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
v = mfC./(dwm.*sC.*(hxhC)*(NC-1));
Re = v.*dwm.*dhC/(mwm);
Pr = pwm;
 

% Colburn and pressure factors
al = sC./(hxhC);
if Re<=1000
j_f = 0.483*(cvl/dhC)^(-0.162)*al^(-0.184)*Re^(-0.536);
f_f = 7.661*(cvl/dhC)^(-0.384)*al^(-0.092)*Re^(-0.712);

elseif Re>=1000
% j_f = 0.242*(1/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% f_f = 1.136*(l/Dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% j_f = 0.242*(cvl/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL

j_f =0.5*(cvl/dhC)^(-0.322)*(tf/dhC)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
f_f = 1.136*(cvl/dhC)^(-0.3)*(tf/dhC)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
end



 h = j_f.*Re.*(Pr.^(1./3)).*kwm./(dhC);
end




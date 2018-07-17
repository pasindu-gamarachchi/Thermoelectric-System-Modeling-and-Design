% Function calculates the convection coefficient along with velocity and 
% Reynolds number for the cold side of the heat exchanger which uses water
% as the working fluid
 function h = coldsideconvcoeff( Tc,mf)
% 
 global hxl hxh hxw N_C tf_C  cvl Tcg

%Test Inputs
% Tc= 25;
% mf = 0.1;
% hxl = 0.16;
% hxh = 0.2;
% hxw = 0.04;
% N_C =5;
% tf_C = 0.0001;
% cvl = 0.02;
% Tcg =5;




Tci= Tc +Tcg;

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

% Thermal Conductivity of Fin Material- Nickel
knf = int(-9.32400932E-11*x^4 + 1.13247863E-07*x^3 +6.33449883E-05*x^2 + -9.47163947E-02*x + 8.13811189E+01);



dw =  symfun(dwf,x);
cw =  symfun(cwf,x);
mw =  symfun(mwf,x);
kw =  symfun(kwf,x);
pw =  symfun(pwf,x);
kn =  symfun(knf,x);

dwm = (dw(Tci) - dw (Tc))./(Tci-Tc);
dwm = double(dwm);

cwm = (cw(Tci) - cw (Tc))./(Tci-Tc);
cwm = double(cwm);

mwm = (mw(Tci) - mw (Tc))./(Tci-Tc);
mwm = double(mwm);

kwm = (kw(Tci) - kw (Tc))./(Tci-Tc);
kwm = double(kwm);

pwm = (pw(Tci) - pw (Tc))./(Tci-Tc);
pwm = double(pwm);

knm = (kn(Tci) - kn (Tc))./(Tci-Tc);
knm = double(knm);

% Intermediate Calculations
s = (hxw - (N_C.*tf_C))./(N_C+1);  % Fin Spacing
dh = (4*s.*(hxh./2))./(s + hxh );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
v = mf./(dwm.*s.*(hxh/2)*(N_C-1))
Re = v.*dwm.*dh/(mwm)
Pr = pwm;


 % Colburn and pressure factors
al = s./(hxh./2);


if Re<=1000
jf = 0.483*(cvl/dh)^(-0.162)*al^(-0.184)*Re^(-0.536);
 f_f = 7.661*(cvl/dh)^(-0.384)*al^(-0.092)*Re^(-0.712);

elseif Re>=1000
% j_f = 0.242*(1/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% f_f = 1.136*(l/Dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% j_f = 0.242*(cvl/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL

jf =0.5*(cvl/dh)^(-0.322)*(tf/dh)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
f_f = 1.136*(cvl/dh)^(-0.3)*(tf/dh)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
end


h = jf.*Re.*(Pr.^(1./3)).*kwm./(dh)

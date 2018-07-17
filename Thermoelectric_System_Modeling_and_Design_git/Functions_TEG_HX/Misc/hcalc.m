% Function calculates the convection coefficient along with Velocity and
% Reynolds Number for the hot side of the heat exchanger which uses air as
% the working fluid 
 function [h v Re] = hcalc(Tcold, mf)

 global hxl hxh hxw N tf   tc sl Thg cvl Tcg

 
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



kb = symfun(knf,x);
% Integral averages

dwm = (dw(Tco) - dw (Tcold))./(Tcg);
dwm = double(dwm);

cwm = (cw(Tco) - cw(Tcold))./(Tcg);
cwm = double(cwm);

mwm = (mw(Tco) - mw (Tcold))./(Tcg);
mwm = double(mwm);

kwm = (kw(Tco) - kw(Tcold))./(Tcg);
ksm = double(kwm);

pwm = (pw(Tco) - pw(Tcold))./(Tcg);
pwm = double(psm);



cam = (ca(Th) - ca (Thi))./(Th-Thi);
cam = double(cam);

mam = (ma(Th) - ma (Thi))./(Th-Thi);
mam = double(mam);

kam = (ka(Th) - ka (Thi))./(Th-Thi);
kam = double(kam);

pam = (pa(Th) - pa (Thi))./(Th-Thi);
pam = double(pam);




% Intermediate Calculations
s = (hxw - (N.*tf))./(N+1);  % Fin Spacing
dh = (4*s.*(hxh./2))./(s + hxh );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
v = mf./(dam.*s.*(hxh/2)*(N-1));
Re = v.*dam.*dh/(mam);
Pr = (cam.*mam)./(kam);
 

% Colburn and pressure factors
al = s./(hxh./2);
if Re<=1000
j_f = 0.483*(cvl/dh)^(-0.162)*al^(-0.184)*Re^(-0.536);
f_f = 7.661*(cvl/dh)^(-0.384)*al^(-0.092)*Re^(-0.712);

elseif Re>=1000
% j_f = 0.242*(1/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% f_f = 1.136*(l/Dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% j_f = 0.242*(cvl/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL

j_f =0.5*(cvl/dh)^(-0.322)*(tf/dh)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
f_f = 1.136*(cvl/dh)^(-0.3)*(tf/dh)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
end



 h = j_f.*Re.*(Pr.^(1./3)).*kam./(dh);
% Re





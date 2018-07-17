% Function calculates the convection coefficient along with Velocity and
% Reynolds Number for the hot side of the heat exchanger which uses water as
% the working fluid , also calculates pressure drop in control volume
 function [h, v, Re, f_f, prD, V_f] = convcoeff_DuctWater(Tc, mf)

 global hxl hxh hxw Tcg 

Tin=Tc;
Tco = Tc + Tcg;

%%  Temperature dependent properties - Water

% Density Water
 syms x % ps cs ms ks ps 
dwf = int(-1.244106E-07*x^4 + 4.224596E-05*x^3 - 7.670379E-03*x^2 + 5.169102E-02*x + 1.000067E+03);

% Specific Heat Capacity of Water
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


dwm = (dw(Tco) - dw (Tc))/(Tco-Tc);
dwm = double(dwm);

cwm = (cw(Tco) - cw(Tc))/(Tco - Tc);
cwm = double(cwm);

mwm = (mw(Tco) - mw(Tc))./(Tco - Tc);
mwm = double(mwm);

kwm = (kw(Tco) - kw(Tc))./(Tco-Tc);
kwm = double(kwm);

pwm = (pw(Tco) - pw(Tc))./(Tco - Tc);
pwm = double(pwm);






% Intermediate Calculations
dh  = (2*hxw*hxl)/(hxl/2 +hxw/2);
v = mf/(dwm*hxw*hxh); % Channel Velocity
% Re2 = v*dam*dh/(mam);
% Re = ((mf/(N-1))*4)/(pi*dh*mam);
Re = v*dwm*dh/(mwm);
Pr = (cwm*mwm)/(kwm);
 

% Colburn and pressure factors
al = hxh/(hxw);

if Re < 2500

    Nu = 7.541*(1- 2.61*al +4.97*al^2 -5.199*al^3 + 2.702*al^4 -0.548*al^5);
else
    
    Nu = 0.026*Re^(0.8)*Pr^(0.3);
    
end

h = Nu*kwm/dh;

if Re < 2100
    f_f = (24/Re)*(1 - 1.3553*al + 1.9467*al^2 - 1.7012*al^3 + 0.9564*al^4 -0.2537*al^5);
    
elseif 2100 <= Re <= 4100
    A = 0.0054;
    B = 2.3*(10^-8);
    m = (-2/3);
    f_f  = A + B*(Re^(-1/m));
    
elseif Re > 4100
    A = 0.00128;
    B  = 0.1143;
    m  = 3.2154;
    f_f  = A + B*(Re^(-1/m));
end
    

prD = 0.5*f_f*hxl*dwm*v^2/(dh);
V_f = mf/dwm;


 end



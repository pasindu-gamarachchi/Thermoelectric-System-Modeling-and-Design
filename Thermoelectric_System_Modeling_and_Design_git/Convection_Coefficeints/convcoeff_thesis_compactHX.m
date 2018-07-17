% Function calculates the convection coefficient along with Velocity and
% Reynolds Number for the hot side of the heat exchanger 
% Pasindu Gamarachchi - Email : pgamarachchi@gmail.com
 function [h, v, Re] = convcoeff(Th, mf)

 global hxl hxh hxw N tf Thg 

Thi = Th - Thg;

%%  Temperature dependent properties

% Density of steam/ Hot side fluid
 syms x % ps cs ms ks ps 
daf = int(2.876602E-12*x^4 - 7.350893E-09*x^3 + 7.284062E-06*x^2 - 3.760334E-03*x + 1.251051E+00);
% Specific Heat Capacity of steam
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03);
% Dynamic viscocity of Steam
maf = int(-1.17230617E-17*x^4 +3.14251436E-14*x^3 -3.87294440E-11*x^2 +4.96283182E-08*x + 1.71301016E-05);
% Thermal Cond of Steam
kaf = int(1.146049E-14*x^4 + 5.825050E-12*x^3 - 4.412201E-08*x^2 + 8.330657E-05*x + 2.396400E-02);  
% Prandtl Number of Steam
paf = int( 1.946038E-13*x^4 - 8.530084E-10*x^3 + 9.958116E-07*x^2 - 3.351574E-04*x + 7.176453E-01);
% Thermal Conductivity of Fin Material
knf = int(-9.32400932E-11*x^4 + 1.13247863E-07*x^3 +6.33449883E-05*x^2 + -9.47163947E-02*x + 8.13811189E+01);

da =  symfun(daf,x);
ca =  symfun(caf,x);
ma =  symfun(maf,x);
ka =  symfun(kaf,x);
pa =  symfun(paf,x);

kb = symfun(knf,x);
% Integral averages


cam = (ca(Th) - ca (Thi))/(Th-Thi);
cam = double(cam);

mam = (ma(Th) - ma (Thi))/(Th-Thi);
mam = double(mam);

kam = (ka(Th) - ka (Thi))/(Th-Thi);
kam = double(kam);

pam = (pa(Th) - pa (Thi))/(Th-Thi);
pam = double(pam);

% Intermediate Calculations
s = (hxw - (N*tf))/(N+1); % Fin Spacing
dh = (4*s*(hxh/2))/(s + hxh );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
v = mf/(dam*s*(hxh/2)*(N-1));
Re = v*dam*dh/(mam);
Pr = (cam*mam)/(kam);
 

% Colburn and pressure factors
al = s/(hxh/2);
if Re<=1000
    j_f = 0.483*(hxl/dh)^(-0.162)*al^(-0.184)*Re^(-0.536);
else 
    j_f =0.242*(hxl/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368);
end

h = j_f*Re*(Pr^(1/3))*kam/(dh);

 end



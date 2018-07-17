% Function calculates the convection coefficient along with Velocity and
% Reynolds Number for the hot side of the heat exchanger which uses air as
% the working fluid , also calculates pressure drop in control volume
 function [h, v, Re, f_f, f_f2, prD, prD_O, prD_3] = convcoeffv3(Th, mf)

 global hxl hxh hxw N tf   tc sl Thg cvl

 Tin=Th;
Thi = Th - Thg;

%%  Temperature dependent properties

% Density of steam/ Hot side fluid
 syms x % ps cs ms ks ps 
dsf = int(4.19200381E-12*x^4 - 1.25128021E-08*x^3 + 1.45079292E-05*x^2 - 8.12253344E-03*x + 2.17634458E+00);
daf = int(2.876602E-12*x^4 - 7.350893E-09*x^3 + 7.284062E-06*x^2 - 3.760334E-03*x + 1.251051E+00);
% Specific Heat Capacity of steam
csf = int( 2.38881298E-08*x^4 - 6.51294808E-05*x^3 + 6.61776463E-02*x^2 - 2.90862268E+01*x + 6.62561202E+03 );
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03);
% Dynamic viscocity of Steam
msf = int(3.60777971E-08*x - 9.93254096E-07);
maf = int(-1.17230617E-17*x^4 +3.14251436E-14*x^3 -3.87294440E-11*x^2 +4.96283182E-08*x + 1.71301016E-05); %int(7.012843E-18*x^4 + 2.116562E-14*x^3 - 3.200137E-11*x^2 + 4.954035E-08*x + 1.717161E-05);
% Thermal Cond of Steam
ksf = int( 1.38151501E-08*x^2 + 6.63499457E-05*x - 2.67362298E-03);
kaf = int(1.146049E-14*x^4 + 5.825050E-12*x^3 - 4.412201E-08*x^2 + 8.330657E-05*x + 2.396400E-02); %int( -5.10944955E-15*x^4 + 2.56059110E-11*x^3 -5.79413473E-08*x^2 + 1.05976047E-04*x -9.82822783E-04);%int( 1.146049E-14*x^4 + 5.825050E-12*x^3 - 4.412201E-08*x^2 + 8.330657E-05*x + 2.396400E-02);
% Prandtl Number of Steam
psf = int( 1.15473535E-11*x^4 - 3.07625489E-08*x^3 + 3.07401737E-05*x^2 - 1.36012594E-02*x + 3.23609080E+00);
paf = int( 1.946038E-13*x^4 - 8.530084E-10*x^3 + 9.958116E-07*x^2 - 3.351574E-04*x + 7.176453E-01);

% Thermal Conductivity of Fin Material
knf = int(-9.32400932E-11*x^4 + 1.13247863E-07*x^3 +6.33449883E-05*x^2 + -9.47163947E-02*x + 8.13811189E+01);

ds =  symfun(dsf,x);
cs =  symfun(csf,x);
ms =  symfun(msf,x);
ks =  symfun(ksf,x);
ps =  symfun(psf,x);

da =  symfun(daf,x);
ca =  symfun(caf,x);
ma =  symfun(maf,x);
ka =  symfun(kaf,x);
pa =  symfun(paf,x);

kb = symfun(knf,x);
% Integral averages

dsm = (ds(Th) - ds (Thi))/(Th-Thi);
dsm = double(dsm);

csm = (cs(Th) - cs (Thi))/(Th-Thi);
csm = double(csm);

msm = (ms(Th) - ms (Thi))/(Th-Thi);
msm = double(msm);

ksm = (ks(Th) - ks (Thi))/(Th-Thi);
ksm = double(ksm);

psm = (ps(Th) - ps (Thi))/(Th-Thi);
psm = double(psm);

dam = (da(Th) - da (Thi))/(Th-Thi);
dam = double(dam);


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
dh = (4*s*(hxh))/(2*s + 2*hxh );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
v = 2*mf/(dam*s*(hxh)*(N-1)); % Channel Velocity
% Re2 = v*dam*dh/(mam);
% Re = ((mf/(N-1))*4)/(pi*dh*mam);
Re = v*dam*dh/(mam);
Pr = (cam*mam)/(kam);
 

% Colburn and pressure factors
al = s/(hxh/2);
% if Re<=2000   %if Re<=1000
%     j_f = 0.483*(hxl/dh)^(-0.162)*al^(-0.184)*Re^(-0.536);
%     f_f = 7.661*(hxl/dh)^(-0.384)*al^(-0.092)*Re^(-0.712);
% % j_f = 0.483*(cvl/dh)^(-0.162)*al^(-0.184)*Re^(-0.536);
% %  f_f2 = 64/Re;
% 
% else 
% % j_f = 0.242*(1/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
%  f_f = 1.136*(hxl/dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% j_f = 0.242*(hxl/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% 
% % j_f =0.5*(cvl/dh)^(-0.322)*(tf/dh)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
% % f_f = 1.136*(cvl/dh)^(-0.3)*(tf/dh)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
% 
% %     j_f =0.5*(hxl/dh)^(-0.322)*(tf/dh)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
% %     f_f = 1.136*(hxl/dh)^(-0.3)*(tf/dh)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
% %     f_f2 = 1.136*(hxl/dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% 
% end
% 
 if Re<=1000   %if Re<=1000
% %     j_f = 0.483*(hxl/dh)^(-0.162)*al^(-0.184)*Re^(-0.536);
    f_f2 = 7.661*(hxl/dh)^(-0.384)*al^(-0.092)*Re^(-0.712);
% % j_f = 0.483*(cvl/dh)^(-0.162)*al^(-0.184)*Re^(-0.536);
% %  f_f2 = 64/Re;
% 
 else 
% % j_f = 0.242*(1/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
f_f2 = 1.136*(hxl/dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% % j_f = 0.242*(hxl/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% 
% % j_f =0.5*(cvl/dh)^(-0.322)*(tf/dh)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
% % f_f = 1.136*(cvl/dh)^(-0.3)*(tf/dh)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
% 
% %     j_f =0.5*(hxl/dh)^(-0.322)*(tf/dh)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
% %     f_f = 1.136*(hxl/dh)^(-0.3)*(tf/dh)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
% %     f_f2 = 1.136*(hxl/dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% 
 end

Nu = 7.541*(1- 2.61*al +4.97*al^2 -5.199*al^3 + 2.702*al^4 -0.548*al^5);


h = Nu*kam/dh;

f_f = (24/Re)*(1 - 1.3553*al + 1.9467*al^2 - 1.7012*al^3 + 0.9564*al^4 -0.2537*al^5);

if f_f > 0.1
    f_f = 0.08;
else
end


if Re > 500
    h = (0.024*Re^(0.8)*(Pr^(0.4)))*kam/dh;
else
end
%  h = j_f*Re*(Pr^(1/3))*kam/(dh);

% prD = (N-1)*2*f_f*cvl*dam*v^2/(dh);
prD = (N-1)*0.5*f_f*cvl*dam*v^2/(dh);
prD_O = dam*v^2*(N-1)*f_f2*cvl/dh;
prD_3 = dam*v^2*(N-1)*f_f*cvl/dh;

 end



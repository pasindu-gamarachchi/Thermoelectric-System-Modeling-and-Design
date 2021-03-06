% Function calculates the heat transfer in a specified control volume for
% the hot side of the heat exchanger
  function h = convcoeff(Th, mf)

 global hxl hxh hxw N tf   tc sl
% % Inputs
% hxh = 0.1; % Height of heat exchanger [m]
% hxw = 0.3; % Width of heat exchanger [m]
% hxl = 5.*(1.8*10^-3); % Length of Control Volume [m]
% F = 5;  % Number of fins
% 
% tf = 0.01; % Fin thickness [m]

%  mfv = 5.*linspace(0.005,0.01,10);
% for i =1:length(mfv)
% mf = mfv(i); % kg/s Mass flowrate
%  Th = 558; % C - Inlet Temperature
% Thi = 495;
 Tin=Th;
Thi = Th - 10;

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
kaf =  int( -5.10944955E-15*x^4 + 2.56059110E-11*x^3 -5.79413473E-08*x^2 + 1.05976047E-04*x -9.82822783E-04);%int( 1.146049E-14*x^4 + 5.825050E-12*x^3 - 4.412201E-08*x^2 + 8.330657E-05*x + 2.396400E-02);
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

dsm = (ds(Th) - ds (Thi))./(Th-Thi);
dsm = double(dsm);

csm = (cs(Th) - cs (Thi))./(Th-Thi);
csm = double(csm);

msm = (ms(Th) - ms (Thi))./(Th-Thi);
msm = double(msm);

ksm = (ks(Th) - ks (Thi))./(Th-Thi);
ksm = double(ksm);

psm = (ps(Th) - ps (Thi))./(Th-Thi);
psm = double(psm);

dam = (da(Th) - da (Thi))./(Th-Thi);
dam = double(dam);

cam = (ca(Th) - ca (Thi))./(Th-Thi);
cam = double(cam);

mam = (ma(Th) - ma (Thi))./(Th-Thi);
mam = double(mam);

kam = (ka(Th) - ka (Thi))./(Th-Thi);
kam = double(kam);

pam = (pa(Th) - pa (Thi))./(Th-Thi);
pam = double(pam);


% Geometries
% he = 0.03 ; % [m] Height of chamber
% delx = 1.8*10^-3; % [m] length of chamber
% w  = (1.8*10^-3)*4; % [m] Width of chamber

% Intermediate Calculations
s = (hxw - (N.*tf))./(N+1);  % Fin Spacing
dh = (2.*s.*(hxh./2))./(s + hxh./2 );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
v = mf./(dam.*s.*hxh*(N));
Re = v.*dam.*dh/(mam);
Pr = (cam.*mam)./(kam);
 

% Colburn and pressure factors
al = s./(hxh./2);
if Re<=1000
j_f = 0.483*(hxl/dh)^(-0.162)*al^(-0.184)*Re^(-0.536);
f_f = 7.661*(hxl/dh)^(-0.384)*al^(-0.092)*Re^(-0.712);

elseif Re>=1000
% j_f = 0.242*(l/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% f_f = 1.136*(l/Dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL

j_f =0.5*(hxl/dh)^(-0.322)*(tf/dh)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
f_f = 1.136*(hxl/dh)^(-0.3)*(tf/dh)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
end

% Convection Coefficent

h = j_f.*Re.*(Pr.^(1./3)).*kam./(dh);

% %  end



% Q = mf.*(Th-Thi).*cam;



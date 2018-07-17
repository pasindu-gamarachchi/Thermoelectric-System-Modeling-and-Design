% Cold Side Heat Exchanger
% Function Calculates Heat Transferred from 
% cold side fins to cold fluid 

% Updated on 04/06
% Working fluid changed to air
 function [Qc Tout h] = HXColdAir(Tcb, Tcf)

global NC tf_C mfC cvl Tcg Tcbg hxlC hxhC hxwC

% Tcb =70;
% Tcf =20;

% All cold side Heat Exchanger Parameters in this file
Tcg  = 8;


Tcbg= 120;




% hxlC = 0.16; % [m]
hxhC = 0.02; % [m]
hxwC = 0.04; % [m]

% PackFraction = (tf_C*NC/hxwC)

%%
% sC = (hxw - (N_C.*tf_C))./(N_C+1); % Fin Spacing

syms x
% 
% cwf = int(5.951927E-06*x^4 - 1.221571E-03*x^3 + 9.793076E-02*x^2 - 3.302581E+00*x + 4.217176E+03); % Specific Heat of Water
% cw = symfun(cwf, x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);


% cwm = (cw(Tcg) - cw(Tcf))./(Tcg-Tcf);
% cwm = double(cwm);
% 
% cp =cwm;

Tb = Tcf + Tcbg;

Pf = 2*(cvl);
Af = cvl*tf_C;


knim = (kni(Tb) - kni(Tb-10))./(10);
knim = double(knim);

%% Start of Cold convection function
Tcold = Tcf;
Tco = Tcold + Tcg;
%%  Temperature dependent properties

 syms x  
 % Density of water/ Hot side fluid
% dwf = int(-1.244106E-07*x^4 + 4.224596E-05*x^3 - 7.670379E-03*x^2 + 5.169102E-02*x + 1.000067E+03);
daf = int(2.876602E-12*x^4 - 7.350893E-09*x^3 + 7.284062E-06*x^2 - 3.760334E-03*x + 1.251051E+00);
% Specific Heat Capacity of water
% cwf = int( 5.951927E-06*x^4 - 1.221571E-03*x^3 + 9.793076E-02*x^2 - 3.302581E+00*x + 4.217176E+03);
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03);
% Dynamic viscocity of Water
% mwf = int(4.916802E-11*x^4 - 1.125020E-08*x^3 + 1.050048E-06*x^2 - 5.463958E-05*x + 1.756519E-03);
maf = int(-1.17230617E-17*x^4 +3.14251436E-14*x^3 -3.87294440E-11*x^2 +4.96283182E-08*x + 1.71301016E-05); 
% Thermal Cond of Water
% kwf = int( 1.934482E-10*x^4 - 4.145917E-08*x^3 - 4.292555E-06*x^2 + 1.765088E-03*x + 5.695880E-01);
kaf =  int( -5.10944955E-15*x^4 + 2.56059110E-11*x^3 -5.79413473E-08*x^2 + 1.05976047E-04*x -9.82822783E-04);
% Prandtl Number of Water
% pwf = int( 4.299125E-07*x^4 - 9.796954E-05*x^3 + 9.010003E-03*x^2 - 4.486901E-01*x + 1.305137E+01);
paf = int( 1.946038E-13*x^4 - 8.530084E-10*x^3 + 9.958116E-07*x^2 - 3.351574E-04*x + 7.176453E-01);


% dw =  symfun(dwf,x);
% cw =  symfun(cwf,x);
% mw =  symfun(mwf,x);
% kw =  symfun(kwf,x);
% pw =  symfun(pwf,x);

da =  symfun(daf,x);
ca =  symfun(caf,x);
ma =  symfun(maf,x);
ka =  symfun(kaf,x);
pa =  symfun(paf,x);



% Integral averages

dam = (da(Tco) - da (Tcold))./(Tcg);
dam = double(dam);

cam = (ca(Tco) - ca(Tcold))./(Tcg);
cam = double(cam);

mam = (ma(Tco) - ma(Tcold))./(Tcg);
mam = double(mam);

kam = (ka(Tco) - ka(Tcold))./(Tcg);
kam = double(kam);

pam = (pa(Tco) - pa(Tcold))./(Tcg);
pam = double(pam);







% Intermediate Calculations
sC = (hxwC - (NC.*tf_C))./(NC+1);  % Fin Spacing
dhC = (4*sC.*(hxhC))./(sC + 2*hxhC );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
v = mfC./(dam.*sC.*(hxhC)*(NC-1));
Re = v.*dam.*dhC/(mam);
Pr = (cam.*mam)./(kam);
 

% Colburn and pressure factors
al = sC./(hxhC);
if Re<=1000
j_f = 0.483*(cvl/dhC)^(-0.162)*al^(-0.184)*Re^(-0.536);
f_f = 7.661*(cvl/dhC)^(-0.384)*al^(-0.092)*Re^(-0.712);

elseif Re>=1000
% j_f = 0.242*(1/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL
% f_f = 1.136*(l/Dh)^(-0.781)*(tf/dh)^(0.534)*Re^(-0.198); % ORIGINAL
% j_f = 0.242*(cvl/dh)^(-0.322)*(tf/dh)^0.089*Re^(-0.368); % ORIGINAL

j_f =0.5*(cvl/dhC)^(-0.322)*(tf_C/dhC)^0.5*Re^(-0.368); % MODIFIED BASED ON FLUENT
f_f = 1.136*(cvl/dhC)^(-0.3)*(tf_C/dhC)^(0.44)*Re^(-0.23); % MODIFIED BASED ON FLUENT
end



 h = j_f.*Re.*(Pr.^(1./3)).*kam./(dhC);





M = (sqrt(h*Pf*knim*Af))*(Tcb -Tcf);
m = sqrt((h*Pf)/(knim*Af));

Qc = M*tanh(m*hxhC)*(NC-1);

Tout = Tcf + Qc/(mfC*cam);
% h
% v
% Re
%Toc = Tcf + (Qc./(mfC*cwm))

% h
% v
% Re
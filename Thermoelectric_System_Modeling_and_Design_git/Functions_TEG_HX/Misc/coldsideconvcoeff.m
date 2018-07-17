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






% Intermediate Calculations
s = (hxw - (N_C.*tf_C))./(N_C+1);  % Fin Spacing
dh = (4*s.*(hxh./2))./(s + hxh );  % Hydraulic Diameter, hxh is divided by 2, because of adiabatic tip
v = mf./(dwm.*s.*(hxh/2)*(N_C-1))
Re = v.*dwm.*dh/(mwm)
Pr = pwm;


 % Colburn and pressure factors
al = s./(hxh./2);


if Re<1000
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

% Natural Convection Cold Side Heat Exchanger
% Function Calculates Heat Transferred from Heat Sink Fins
% Calculations based on source : [1] Joo, Y. and Kim, S., 2015, "Comparison of thermal performance between plate-fin and pin-fin heat sinks in natural convection", International Journal of Heat and Mass Transfer, 83, pp. 345-356.
% Heat transfer from Fins, also calculates the Tip Temperature


   function [QT, Q_b, h_f, h_b, TEMC, TTip, et_f, Sh, Sv,k, h_f1,h_f2,h_f3,h_f4] = HXFreePinFin(Tcb, Tamb)


 global NCx NCy  cvl  hxlC  hxwC  Cb  fd fh finmat


g = 9.81;
Boltz = 5.67*10^-8;
Err = (10^-4);


syms x

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

kAlpf = int( 7.00000000E-09*x^4 - 6.18933333E-06*x^3 + 1.25071800E-03*x^2 + 1.47795093E-02*x + 2.28807284E+02);
kAlp = symfun(kAlpf,x);

kCupf = int(5.25E-15*x^6 - 1.34E-11*x^5 + 1.26E-08*x^4 - 5.19E-06*x^3 + 8.71E-04*x^2 - 9.61E-02*x + 3.99E+02);
kCup = symfun(kCupf,x);

ep_Cu =0.28;

%  Temperature dependent properties

 syms x  
 % Density of water/ Hot side fluid

 daf = int(2.876602E-12*x^4 - 7.350893E-09*x^3 + 7.284062E-06*x^2 - 3.760334E-03*x + 1.251051E+00);
% Specific Heat Capacity of water
%  cwf = int( 5.951927E-06*x^4 - 1.221571E-03*x^3 + 9.793076E-02*x^2 - 3.302581E+00*x + 4.217176E+03);
 caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03);
% Dynamic viscocity of Water
%  mwf = int(4.916802E-11*x^4 - 1.125020E-08*x^3 + 1.050048E-06*x^2 - 5.463958E-05*x + 1.756519E-03);
 maf = int(-7.012843E-18*x^4 + 2.116562E-14*x^3 - 3.200137E-11*x^2 + 4.954035E-08*x + 1.717161E-05);
%  (-1.17230617E-17*x^4 +3.14251436E-14*x^3 -3.87294440E-11*x^2 +4.96283182E-08*x + 1.71301016E-05); 
% Thermal Cond of Water
%  kwf = int( 1.934482E-10*x^4 - 4.145917E-08*x^3 - 4.292555E-06*x^2 + 1.765088E-03*x + 5.695880E-01);
kaf =  int( 1.146049E-14*x^4 + 5.825050E-12*x^3 - 4.412201E-08*x^2 + 8.330657E-05*x + 2.396400E-02);
% Prandtl Number of Water
%  pwf = int( 4.299125E-07*x^4 - 9.796954E-05*x^3 + 9.010003E-03*x^2 - 4.486901E-01*x + 1.305137E+01);
 paf = int( 1.946038E-13*x^4 - 8.530084E-10*x^3 + 9.958116E-07*x^2 - 3.351574E-04*x + 7.176453E-01);
% Thermal Diffusivity of Air
faf = int( 1.093417E-16*x^4 - 1.924826E-13*x^3 + 1.754598E-10*x^2 + 1.375914E-07*x + 1.847226E-05);

% Kinematic Viscocity of Air
vaf = int( 9.471589E-18*x^4 - 3.915231E-14*x^3 + 1.054041E-10*x^2 + 8.961175E-08*x + 1.340998E-05);


da =  symfun(daf,x);
ca =  symfun(caf,x);
ma =  symfun(maf,x);
ka =  symfun(kaf,x);
pa =  symfun(paf,x);
fa =  symfun(faf,x); 
va =  symfun(vaf,x);

kAlpm = (kAlp(Tcb) - kAlp(Tcb-0.1))./(0.1);
kAlpm = double(kAlpm);

kCupm = (kCup(Tcb) - kCup(Tcb-0.1))./(0.1);
kCupm = double(kCupm);



if finmat == 'Cu'
    k = kCupm;
    
elseif finmat == 'Al'
    
    k = kAlpm;
else
end

    








% Start of Cold convection function



Tfilm = (Tamb + Tcb)/2;

% Air
dam = (da(Tfilm) - da (Tfilm-0.01))./(0.01);
dam = double(dam);

cam = (ca(Tfilm) - ca(Tfilm-0.01))./(0.01);
cam = double(cam);

mam = (ma(Tfilm) - ma(Tfilm-0.01))./(0.01);
mam = double(mam);

kam = (ka(Tfilm) - ka(Tfilm-0.01))./(0.01);
kam = double(kam);

pam = (pa(Tfilm) - pa(Tfilm-0.01))./(0.01);
pam = double(pam);

fam = (fa(Tfilm) - fa(Tfilm-0.01))./(0.01);
fam = double(fam);

vam = (va(Tfilm) - va(Tfilm-0.01))./(0.01);
vam = double(vam);

TfilmK = Tfilm +273.15;
Be = 1/TfilmK;



% Intermediate Calculations
sCx = (hxwC - (NCx*fd))/(NCx-1); % Fin Spacing, actual spacing between pins x-dir
sCy = (hxlC - (NCy*fd))/(NCy-1); % Fin Spacing, actual spacing between pins y-dir


if sCx <= 0 || sCy <=0
    fprintf('Spacing less than 0\n')
    NCx
    sCx
    NCy
    sCy
    return
end

PackFrac = ((pi/4)*(fd^2))*(NCx*NCy)/(hxwC*hxlC);

A_Base = (hxwC*hxlC) -((pi/4)*(fd^2)*NCx*NCy);

% Base Heat Transfer Calculations
Ra_b = (Tcb-Tamb)*g*Be*(pam)*(hxlC^3)/(vam^2); % Raleigh Number for Base Area
Nu_b = 0.59*(Ra_b^(0.25)); % Nusselt Number for Base Area
h_b = (Nu_b*kam)/(hxlC); % Convection Coeff from Base Area 
Q_b = A_Base*(h_b)*(Tcb-Tamb); % Heat transfer from Base Area

% El = g*Be*(Tcb-Tamb)*(pam)*(sC^4)/(hxlC*vam^2) % Elenbaas Number


et_fin = 0.01; % Fin Efficiency

% Sh_a = (hxwC - fd)/(NCx -1);
% Sh = Sh_a/2;

Sv = (hxlC - fd)/(NCy -1);
Sh = (hxlC*hxwC)/(Sv*NCx*NCy);
 
Perm  = (4*Sh*Sv - pi*fd^2)/48;
Sh_s = Sh/(fd^(3/4) * hxlC^(1/4));
Gr_L  = (g*Be*et_fin)*(Tcb -Tamb)*(hxlC^3)/(vam^2);
Ra_f4 =  g*Be*et_fin*(Tcb - Tamb)*(fd^3)/(fam*vam);

h_f1 =(Sh*Sv/(pi*fd*hxlC))*((Perm/vam)* dam*cam*g*Be*et_fin*(Tcb - Tamb)); %((hxwC*Perm)/(pi*fd*NCx*NCy*mam))*(dam^2 *cam *g*Be)*et_fin*(Tcb - Tamb);
h_f2 = (kam/hxlC)*((0.3669*Sv/fd -0.0494))*(Gr_L^(1/4));
h_f3 = (kam/fd)*( 2.132*Sh_s - 0.4064)*(g*Be*et_fin*(Tcb - Tamb)*(1/(fd*vam*fam)))^(0.188);
h_f4 = (kam/fd)*(0.85*Ra_f4^(0.188));
h_f = (((h_f1^(-1.3) + h_f2^(-1.3) + h_f3^(-1.3))^(-8/-1.3)) + h_f4^-8)^(-1/8);

m = sqrt(4*h_f/(k*fd));
et_f = tanh(m*fh)/(m*fh);




%%

while abs(et_fin -et_f) >  Err
    et_fin = et_fin +(Err*0.8);
    
    Gr_L  = (g*Be*et_fin)*(Tcb -Tamb)*(hxlC^3)/(vam^2);
    Ra_f4 =  g*Be*et_fin*(Tcb - Tamb)*(fd^3)/(fam*vam);

    h_f1 = (Sh*Sv/(pi*fd*hxlC))*((Perm/vam)* dam*cam*g*Be*et_fin*(Tcb - Tamb)); %((hxwC*Perm)/(pi*fd*NCx*NCy*mam))*(dam^2 *cam *g*Be)*et_fin*(Tcb - Tamb);
    h_f2 = (kam/hxlC)*((0.3669*Sv/fd -0.0494))*(Gr_L^(1/4));
    h_f3 = (kam/fd)*( 2.132*Sh_s - 0.4064)*(g*Be*et_fin*(Tcb - Tamb)*(1/(fd*vam*fam)))^(0.188);
    h_f4 = (kam/fd)*(0.85*Ra_f4^(0.188));
    h_f = (((h_f1^(-1.3) + h_f2^(-1.3) + h_f3^(-1.3))^(-8/-1.3)) + h_f4^-8)^(-1/8);

    m = sqrt(4*h_f/(k*fd));
    et_f = tanh(m*fh)/(m*fh);
    
    if et_fin >1
        fprintf('Error!\n')
        return
    end

end




Q_f = h_f*et_f*(pi*fd*fh)*(NCx*NCy)*(Tcb- Tamb);

QT = Q_f + Q_b;
% % Nu =  (576/((Ra)^2) + 2.873/((Ra)^0.5))^(-0.5)%(1/24)*Ra*(al)*(1 - exp(-35/(Ra*(al))))^(0.75); % 
% 
TTip = (Tcb-Tamb)/(cosh(m*fh)) + Tamb;

Rb = Cb/(k*cvl*hxwC);
% 
TEMC = Rb*QT + Tcb;

Nu_Spar = QT/(hxlC*(Tcb-Tamb)*kam);
Ra_Spar = g*Be*(Tcb -Tamb)*(hxlC^3)/(fam*vam);

 end

% Q_f
% Q_b
% QT
% h_f
% TEMC
% TTip
% et_fin
% PackFracPer = PackFrac*100
% Sh_mm = Sh*1000
% Sv_mm = Sv*1000
% 
% toc

% 
% et_finv = 0.01:0.01:1;
% 
% for i = 1: length(et_finv)
%     
%     Gr_L  = (g*Be*et_finv(i))*(Tcb -Tamb)*(hxlC^3)/(vam^2);
%     Ra_f4 =  g*Be*et_finv(i)*(Tcb - Tamb)*(fd^3)/(fam*vam);
%     h_f1 = (Sh*Sv/(pi*fd*hxlC))*((Perm/vam)* dam*cam*g*Be*et_finv(i)*(Tcb - Tamb)) 
%     h_f2 = (kam/hxlC)*((0.3669*Sv/fd -0.0494))*(Gr_L^(1/4));
%     h_f3 = (kam/fd)*( 2.132*Sh_s - 0.4064)*(g*Be*et_finv(i)*(Tcb - Tamb)*(1/fd*vam*fam))^(0.188);
%     h_f4 = (kam/fd)*(0.85*Ra_f4^(0.188));
%     h_f = (((h_f1^(-1.3) + h_f2^(-1.3) + h_f^(-1.3))^(-8/-1.3)) + h_f4^-8)^(-1/8);
%     
%     h_fv(i) = h_f;
%     h_f1v(i) = h_f1;
%     h_f2v(i) = h_f2;
%     h_f3v(i) = h_f3;
%     h_f4v(i) = h_f4;
% 
% end
% 
% et_finv(end)
% h_fv(end)
% close all
% figure (1)
% plot(et_finv, h_fv, 'b.')
% hold on
% plot(et_finv, h_f1v, 'k.')
% plot(et_finv, h_f2v, 'r.')
% plot(et_finv, h_f3v, 'g.')
% plot(et_finv, h_f4v, 'm.')
% 



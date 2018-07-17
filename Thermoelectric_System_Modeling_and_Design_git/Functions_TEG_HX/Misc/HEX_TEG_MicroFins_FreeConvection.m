% 
clc
clear 
close all
format long
tic

% Version Developed 05/23/2016
% Cold Side Heat Sink Natural Convection, uses micro fin heat sink
% TEM Hot Side set to a constant Temperature of 37C
% Separate Function Used for Cold Side Heat Exchanger
% 06/01 spacing equation updated



global hxl  TempC cvl fd  fh NCy  NCx Tcg Tcbg hxhC hxwC hxlC Cb hxw


Tamb= 20;


% Cold Side Dimensions

% tf_C = 0.01*(10^-3); 

% hxhC = 0.01;  %m
hxwC = 0.04; %m
hxl  = 0.04; %m
hxlC = hxl;
Cb = 0.5*10^-2; % m Cold Base Thickness

% TC Dimensions
tc = (1.8*10^-3); % Leg width
sl = 0.4*10^-3; % Spacing between unicouple legs
cs = 1*10^-3; % Spacing between Couples
cw = tc*2 +sl;
cvl = cw +cs;

hxw = hxwC;
cvl = hxl; % Using TEBlock
cv  = floor(hxl./(cvl)); % Number of control volumes

% Fin dimensions and Parameters

fd = 10*(10^-6);
fh = 5*(10^-3);
NCx = 200;
NCy = 200;

SA= pi*fd*fh*NCx*NCy % Total Fin Surface Area

FinA = (NCx*NCy)*(pi/4)*(fd^2) % Fin Cross-Sectional Area

PackFrac = FinA/(hxwC*hxl) % Pack Fraction based on cross-sectional area of fins and total base area

sC = (hxwC -(fd*NCx))/(NCx -1)
pinfin_cond = sC/(68*10^-9)

Rat = sC/fh

AE = SA/(hxwC*hxl)

Th = 37; % Hot Side Temperature

% Guesses and Error
Err = 0.05; % Convergence Error

Tcg = 0.004; % [C] Guess for temperature of TEMC
Tcbg = 7.5; %  [C] Guess for temperature of Fin Base
ErrMult = 8.5;



syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);



% Input Checks

% PackFracColdSide = (tf_C*NC/hxwC)
% sC = (hxwC - (NC.*tf_C))./(NC-1)
Tb = Tamb + Tcbg;
TempC = Tb+0.04; % [C]

QTE = TEBlock(Th, TempC);

% Tb =40;

[Qc hC] = HXMicrofins(Tb, Tamb);

Rth = (Tb-Tamb)/Qc;
    %%


    while ( abs(Qc -QTE)> Err  )
        

         QTE = TEBlock(Th, TempC);
        [Qc hC TEMC ] = HXMicrofins(Tb, Tamb);%mean([Tcf(j), Tcf(j+1)]) );
       
        TempC = TEMC;
        if Qc -QTE> Err & Qc -QTE< (Err*ErrMult)   %Qc -QTE> Err & j ==1
          Tb = Tb -(Err*0.05);
        elseif Qc - QTE> (Err*ErrMult) 
            Tb = Tb -(Err*8.0);
        elseif (Qc -QTE) < (Err*-1) & Qc -QTE > (Err*-ErrMult)  %abs(Qc -QTE) > (Err*1) 
            Tb = Tb +(Err*0.05);
        elseif (Qc -QTE) < (Err*-ErrMult)
            Tb = Tb +(Err*8.0);
        end
        
%         
         
        diffcs  = Qc-QTE
        Tb
%         
%         progress = (j-1)/cv
    end
 
 %%
%  Sp_Cond = Ra*Rat;
 Rth = (Tb-Tamb)/Qc
 
%  fprintf('Thermal Resistance is %
 fprintf('Thermal Resistance is %d C/W\n', Rth)
fprintf('TEM Cold Side Temp is %d C\n', TEMC)
fprintf('Cold Side Fin Base Temp is %d C\n', Tb)
fprintf('Heat Transferred is %d W\n', Qc)
fprintf('Convection Coeff is %d W/m^2-C\n', hC)
fprintf('Fin Spacing is %d mm\n', sC)
fprintf('Fin Surface Area is %d m^2\n', SA)
fprintf('Fin Cross-Sectional Area is %d m^2\n', FinA)
fprintf('Packing Fraction is %d\n', PackFrac)
fprintf('Fin Spacing is %d m\n', sC)
fprintf('Spacing Condition is %d\n', Sp_Cond)
fprintf('Area Enhancement is %d\n', AE)
pinfin_cond = sC/(68*10^-9)

    




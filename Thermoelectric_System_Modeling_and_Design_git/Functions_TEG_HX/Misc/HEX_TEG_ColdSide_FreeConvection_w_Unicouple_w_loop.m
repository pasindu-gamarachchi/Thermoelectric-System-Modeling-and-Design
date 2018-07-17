% 
clc
clear 
close all
format long
tic

% Version Developed 05/10/2016
% Cold Side Heat Sink Natural Convection
% TEM Hot Side set to a constant Temperature of 37C
% Separate Function Used for Cold Side Heat Exchanger
% Modified 05/18/2016 TE Block replaced by Unicouple
% Inital  HX_TEG_Coldside_FreeConvection model modified to include loop for
% optimization.


global hxl  TempC cvl NC tf_C Tcg Tcbg hxhC hxwC hxlC Cb hxw hn ln tn hp lp tp      

Th = 300; % Hot Side Temperature
Tamb= 20;

        
%%


% Cold Side Dimensions
NC = 7;
tf_C = 1.4*(10^-3);


hxhC = 0.12;  %m
hxwC = 0.04; %m
hxl  = 0.04; %m
hxlC = hxl;
Cb = 3*10^-2; % m Cold Base Thickness

% TC Dimensions

hn = 1.69*10^-3;
% hn = 0.2*(10^-3);
ln = 1.3*(10^-3);
tn = 1.3*(10^-3);
hp = hn;
lp = 1.3*(10^-3);
tp = 1.3*(10^-3);

tc = tn; % Leg width
sl = 0.1*10^-3; % Spacing between unicouple legs
cs = 0.1*10^-3; % Spacing between Couples
cw = tc*2 +sl;
cvl = cw +cs;

F = 13; % Couples Per Control Volume
E = floor(hxwC./(cw + cs));

if F>E
    fprintf('Number of Couples per Contol Volume too large\n')
    return
else
end

G = 4;
H  = floor(hxl./(cvl)); % Number of unicouples along fin direction

if G>H
    fprintf('Number of Couples per Contol Volume too large\n')
    return
else
end


hxw = hxwC;
cvl = hxl; 

Th = 300; % Hot Side Temperature

% Guesses and Error
Err = 0.001; % Convergence Error

Tcg = 3; % [C] Guess for temperature of TEMC
Tcbg = 137; %  [C] Guess for temperature of Fin Base
ErrMult = 5.5;
DiffLim = 500; % If the convergence difference is larger than this number*Err, change Temp by Err Mult 
ErrMult = 500;
TempCh = 0.5*Err;

if TempCh > Err
    fprintf('Unlikely to converge\n')
    return
end



syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);


FillFactor = (tn*ln*2*G*F)/(hxwC*hxlC)

% Input Checks

PackFracColdSide = (tf_C*NC/hxwC)
sC = (hxwC - (NC.*tf_C))./(NC-1)
% sCm(i,j) = sC;
Tb = Tamb + Tcbg;
TempC = Tb+3; % [C]

[QHT, QCT, P, n, Volt, I, Qm] = TEModule(Tb,TempC);

QTE = QCT
% QTE = TEBlock(Th, TempC);

% Tb =40;

[Qc, Q_b, h_f, h_b ,TEMC ,TTip ,et_f, Qr] = HXFreeConvV4(Tb, Tamb);

% Rth = (Tb-Tamb)/Qc
    %%


    while ( abs(Qc -QTE)> Err  )
        

         [QHT, QCT, P, n, Volt, I, Qm] = TEModule(Tb,TempC);
         QTE = Q*F*G;
        [Qc, Q_b, h_f, h_b ,TEMC ,TTip ,et_f, Qr] = HXFreeConvV4(Tb, Tamb);%mean([Tcf(j), Tcf(j+1)]) );
       
        TempC = TEMC;
        if Qc -QTE> Err & Qc -QTE< (Err*DiffLim)   %Qc -QTE> Err & j ==1
          Tb = Tb -(TempCh);
        elseif Qc - QTE> (Err*DiffLim) 
            Tb = Tb -(ErrMult*TempCh);
        elseif (Qc -QTE) < (Err*-1) & Qc -QTE > (Err*-DiffLim)  %abs(Qc -QTE) > (Err*1) 
            Tb = Tb +(TempCh);
        elseif (Qc -QTE) < (Err*-ErrMult)
            Tb = Tb +(ErrMult*TempCh);
        end
        
%         
         
        diffcs  = Qc-QTE
        Tb;
%         
%         progress = (j-1)/cv
    end
 
 %%
 Rth = (TEMC-Tamb)/(mean([Qc,QTE]))
Rth_TEM = (Th-TEMC)/(mean([Qc,QTE]))
 TotalP = P*G*F;
 Eff = (TotalP*100)/QTE;
 


toc

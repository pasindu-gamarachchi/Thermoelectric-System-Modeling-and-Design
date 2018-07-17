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


global hxl cvl NC tf_C Tcg Tcbg hxhC hxwC hxlC Cb  hn ln tn hp lp tp      

Th =200; % Hot Side Temperature
Tamb= 20;

hnv =  0.4*(10^-3): (0.05*(10^-3)): 2.0*(10^-3);
leng = length(hnv)
%%


%%HX Dimensions
% Cold Side Dimensions
NC = 7;
tf_C = 1.50*(10^-3);

hxhC = 0.15;  %m
hxwC = 0.04; %m
hxl  = 0.04; %m
hxlC = hxl;
Cb = 1*10^-2; % m Cold Base Thickness
cvl = hxl; 
PackFracColdSide = (tf_C*NC/hxwC)
sC = (hxwC - (NC.*tf_C))./(NC-1)


% TC Dimensions

% hn = 2.4*10^-3;
% hn = 0.2*(10^-3);
ln = 1.8*(10^-3);
tn = 1.8*(10^-3);
% hp = hn;
lp = 1.8*(10^-3);
tp = 1.8*(10^-3);


tc = tn; % Leg width
sl = 0.1*10^-3; % Spacing between unicouple legs
cs = 0.1*10^-3; % Spacing between Couples
cw = tc*2 +sl;
cvl2 = cw +cs;

F = 7; % Couples Per Control Volume
E = floor(hxwC./(cw + cs));

if F>E
    fprintf('Number of Couples per Contol Volume too large\n')
    return
else
end

% Gv =  6:1:12;
G = 7;
H  = floor(hxl./(cvl2)); % Number of unicouples along fin direction
if G>H
    fprintf('Number of Couples per Contol Volume too large\n')
    return
else
end

C = F*G;
FillFactor = (tn*ln*2*C)/(hxwC*hxlC)

% Guesses and Error
Err = 10^-4; % Convergence Error
Tcg = 1; % [C] Guess for temperature of TEMC
Tcbg =155; %  [C] Guess for temperature of Fin Base
ErrMult = 3.5;
DiffLim = 800; % If the convergence difference is larger than this number*Err, change Temp by Err Mult 
ErrMult = 500;
TempCh = 0.3*Err;

if TempCh > Err
    fprintf('Unlikely to converge\n')
    return
end



% syms x
% caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
% ca =  symfun(caf,x);
% 
% knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
% kni =  symfun(knif,x);

%%

for i = 1:length(hnv)

    hn = hnv(i);
    hp = hn;


    if i ==1
        Tb = Tamb + Tcbg;
        TempC = Tb + Tcg;
    else
    end

    [QHT, QCT, P, n, Volt, I, Qm] = TEModuleDOE_HH(Th,TempC,C);

    QTE = QCT;
% QTE = TEBlock(Th, TempC);

% Tb =40;

    [Qc, Q_b, h_f, h_b ,TEMC ,TTip ,et_f, Qr] = HXFreeConvV4(Tb, Tamb);

    %%

    while ( abs(Qc -QTE)> Err  )
        
        [QHT, QCT, P, n, Volt, I, Qm] = TEModuleDOE_HH(Th,TempC,C);
        QTE = QCT;
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
                 
        diffcs  = Qc-QTE
        Tb;
%         
        
    end
 
    progress = i/length(hnv)
 %%
    RthHS = (TEMC-Tamb)/(mean([Qc,QTE]))
    Rth_TEM = (Th-TEMC)/(mean([QHT,QCT]))
    TotalP = P;
    Eff = (TotalP*100)/QTE;

    Qcv(i) =Qc;
    Q_bv(i) = Q_b;
    h_fv(i) = h_f;
    h_bv(i) = h_b;
    TTipv(i) = TTip;
    et_fv(i) = et_f;
    Qrv(i) = Qr;

    TotalPv(i) = P;
    RthHSv(i) = RthHS;
    Rth_TEMv(i) = Rth_TEM;
    QHTv(i) = QHT;
    QCTv(i) = QCT;
    nv(i) =n;
    Voltv(i) = Volt;
    Iv(i) = I;
    Qmv(i) = Qm;
    Tbv(i) =Tb;
    TEMCv(i) = TEMC;

 
end

Ratio = Rth_TEMv./RthHSv;
toc
num2str(Th)
filename = ['PD_DOE_HH_FlatPlate', num2str(Th),'.mat'];
save(filename)


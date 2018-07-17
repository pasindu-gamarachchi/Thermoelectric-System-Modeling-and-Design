% 
clc
clear all
close all
format long
tic

% Version Developed 06/07/2016
% Cold Side Heat Sink Micro-finned, with numerous microwires
% TEM Hot Side set to a constant Temperature of 37C
% Separate Function Used for Cold Side Heat Exchanger
% Modified 05/18/2016 TE Block replaced by Unicouple
% Updated 06/16
% Updated 07/13 Thermal Resistance of human skin added to model

tic
global   cvl Tcg hxwC hxlC Cb NCx NCy fd fh  hn ln tn hp lp tp                              %  hxl  TempC cvl NC tf_C Tcg Tcbg hxhC hxwC hxlC Cb hxw

Th = 37; % Hot Side Temperature
Tamb= 22;


hnv = (0.6*10^-3):(0.020*10^-3):1.28*10^-3
length(hnv)

%%
% Guesses and Error
Err = (10^-6);% 0.01; % Convergence Error
Tcg = 0.04; % [C] Guess for temperature of TEMC
Tcbg = 5.5; %  [C] Guess for temperature of Fin Base
DiffLim = 100; % If the convergence difference is larger than this number*Err, change Temp by Err Mult 
ErrMult = 200;
TempCh = 0.9*Err;
Thg = 3.7;

if TempCh > Err
    fprintf('Unlikely to converge\n')
    return
end



%% HX Dimensions
hxwC = 0.04; %m
hxl  = 0.04; %m
hxlC = hxl;
Cb = 5*10^-3; % m Cold Base Thickness
NCx  = 16;
NCy  = 2; 
fd = 1.0*(10^-3);
fh  = 10*(10^-3);


HXPackFrac = (100)*(pi/4)*(fd^2)*(NCx*NCy)/(hxwC*hxl);

% Leg Dimensions

% hn = 0.2*(10^-3);
ln = 0.6*(10^-3);
tn = 0.6*(10^-3);

lp = 0.6*(10^-3);
tp = 0.6*(10^-3);

% TC Dimensions
tc = tn; % Leg width
sl = 0.4*10^-3; % Spacing between unicouple legs
cs = 2*10^-3; % Spacing between Couples
cw = tc*2 +sl;
tcw = cw +cs;

F = 6; % Couples Per Control Volume
E = floor(hxwC./(tcw));

if F>E
    fprintf('Number of Couples per Contol Volume too large\n')
    return
else
end

G = 3;
H  = floor(hxl./(tcw)); % Number of unicouples along fin direction

if G>H
    fprintf('Number of Couples per Contol Volume too large\n')
    return
else
end

C = G*F

hxw = hxwC;
cvl = hxl; % Using TEBlock






  for a= 1:length(hnv)

hn = hnv(a);
hp = hn;
% Input Checks

% PackFracColdSide = (tf_C*NC/hxwC)
% sC = (hxwC - (NC.*tf_C))./(NC-1)
if a ==1
     Tb = Tamb + Tcbg;
     TEMH = Th - Thg;
else
end

TempC = Tb+ Tcg; % [C]


QHand = ContactRes(Th, TEMH);
[QTEH, QTEC, P, V, I, Qmean ] = TEModule(TEMH,TempC, C);


% QTE = TEBlock(Th, TempC);

% Tb =40;


 [Qc, Q_b, h_f, h_b, TEMC, TTip, et_f, Sh, Sv] = HXFreePinFin(Tb, Tamb);

LegPackFrac = (ln*tn*2*C)/(hxlC*hxwC)
% Rth = (Tb-Tamb)/Qc;
  
%%

i = 1;
    while  abs(Qc -QTEC)> Err  ||  abs(QHand - QTEH) > Err
        

         [QTEH, QTEC, P, N, V, I, Qmean] = TEModule(TEMH, TempC, C);
     
        [Qc, Q_b, h_f, h_b, TEMC, TTip, et_f, Sh, Sv] = HXFreePinFin(Tb, Tamb);

        
        TempC = TEMC;
        QHand = ContactRes(Th, TEMH);

        
        if Qc -QTEC> Err & Qc -QTEC< (Err*DiffLim)   %Qc -QTE> Err & j ==1
          Tb = Tb -(TempCh);
        elseif Qc - QTEC> (Err*DiffLim) 
            Tb = Tb -(ErrMult*TempCh);
        elseif (Qc -QTEC) < (Err*-1) & Qc -QTEC > (Err*- DiffLim)  %abs(Qc -QTE) > (Err*1) 
            Tb = Tb +(TempCh);
        elseif (Qc -QTEC) < (Err*- DiffLim)
            Tb = Tb +(ErrMult*TempCh);
        end
        
        if QHand - QTEH> Err & QHand -QTEH< (Err*DiffLim)   %Qc -QTE> Err & j ==1
            TEMH = TEMH +(TempCh);
        elseif QHand - QTEH> (Err*DiffLim) 
            TEMH  = TEMH  + (ErrMult*TempCh);
        elseif (QHand -QTEH ) < (Err*-1) & QHand -QTEH > (Err*- DiffLim)  %abs(Qc -QTE) > (Err*1) 
            TEMH = TEMH -(TempCh);
        elseif (QHand - QTEH) < (Err*- DiffLim)
            TEMH = TEMH -(ErrMult*TempCh);
        end
        
        if TEMH < TEMC
            fprintf('Error!\n')
            return
        else
        end
        
    

        
        diffhs  = QHand - QTEH
        diffcs  = Qc-QTEC
        
        diffhs_prev = QHand;
        diffcs_prev = Qc;

%         if diffhs <diffhs_prev && diffcs < diffcs_prev
%             fprintf('converging\n')
%             progress = (a-1)*100/length(hnv)
%             
%         else
%             frprintf('Possible Error\n')
%         end
%         
%         diffhs_prev = diffhs;
%         diffcs_prev = diffcs;
%         i = i +1;
     

    end
    
 progress = a/length(hnv)
    
 RthHS = (TEMC-Tamb)/mean([Qc,QTEC]);
 RthHand = (Th-TEMH)/mean([QHand,QTEH]);
 RthTEM = (TEMH -TEMC)/Qmean;
 TotalP = P;
 Eff = (TotalP*100)/QTEH
%  Rthe = RthHSv + RthHandv;
% Ratio = Rth_TEMv./RthHe;
 
 Iv(a) =I;
 Vv(a) =V;
 RthHSv(a) =RthHS;
 Rth_TEMv(a) =RthTEM;
 RthHandv(a) = RthHand;
 TotalPv(a) = TotalP;
 Effv(a) =Eff;
 Qmeanv(a) = Qmean;
 TEMHv(a) = TEMH;
 TEMCv(a) = TEMC;
 QTEHv(a) = QTEH;
 QTECv(a) = QTEC;
%  Grv(a) =Gr;
 
 
 end
%%
Rthe = RthHSv + RthHandv;
Ratio = Rth_TEMv./Rthe;

[Y, In] = max(TotalPv)
Rat_MaxPower =Ratio(In)
HeatIn_MaxPower = QTEHv(In)
LegH_MaxPower = hnv(In)
Volt_MaxPower = Vv(In)
HotT_MaxPower = TEMHv(In)
ColdT_MaxPower = TEMCv(In)

%  Rthv(a) =Rth;
%  Rth_TEMv(a) = Rth_TEM;
%  TotalPv(a) = TotalP;
%  Effv(a) = Eff;
%  TEMCv(a) = TEMC;
%  Qcv(a) =Qc;
%  hcv(a) =hC;
% %  LHV(a) = LH;
% %  LCV(a) = LC;
% %  

%   
%  end
 
 toc
% 
%  Per = '%';
%  fprintf('Thermal Resistance is %
%  fprintf('Thermal Resistance of Heat Sink is %d C/W\n', Rth)
% fprintf('TEM Cold Side Temp is %d C\n', TEMC)
% fprintf('Cold Side Fin Base Temp is %d C\n', Tb)
% fprintf('Heat Transferred is %d W\n', Qc)
% fprintf('Convection Coeff is %d W/m^2-C\n', hC)
% % fprintf('Fin Spacing is %d mm\n', sC)
% fprintf('Total Power is %d W\n', TotalP)
% fprintf('Efficiency is %d %s \n', Eff, Per)
%  fprintf('Thermal Resistance of TEM is %d C/W\n', Rth_TEM)
%  
% end




% syms x
% caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
% ca =  symfun(caf,x);
% 
% knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
% kni =  symfun(knif,x);


% 
clc
clear 
close all
format long
tic

% Version Developed 06/07/2016
% Cold Side Heat Sink Micro-finned, with numerous microwires
% TEM Hot Side set to a constant Temperature of 37C
% Separate Function Used for Cold Side Heat Exchanger
% Modified 05/18/2016 TE Block replaced by Unicouple
% Updated 06/16

tic
global   cvl Tcg hxwC hxlC Cb NCx NCy fd fh  hn ln tn hp lp tp Err TempCh                              %  hxl  TempC cvl NC tf_C Tcg Tcbg hxhC hxwC hxlC Cb hxw

Th = 37; % Hot Side Temperature
Tamb= 22;


hnv = (1.3*10^-3):(0.1*10^-3):2.3*10^-3;
%%


% Guesses and Error
Err = 0.001;% 0.01; % Convergence Error
Tcg = 0.01; % [C] Guess for temperature of TEMC
Tcbg = 8.2; %  [C] Guess for temperature of Fin Base
DiffLim = 100; % If the convergence difference is larger than this number*Err, change Temp by Err Mult 
ErrMult = 200;
TempCh = 0.9*Err;


if TempCh > Err
    fprintf('Unlikely to converge\n')
    return
end


syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

hxwC = 0.04; %m
hxl  = 0.04; %m
hxlC = hxl;
Cb = 5*10^-3; % m Cold Base Thickness
NCx = 13;
NCy = 13;
fd = 69*(10^-6);
fh = 5*(10^-3);





 for a= 1:length(hnv)






% Leg Dimensions
hn = hnv(a);
% hn = 0.2*(10^-3);
ln = 1.8*(10^-3);
tn = 1.8*(10^-3);
hp = hn;
lp = 1.8*(10^-3);
tp = 1.8*(10^-3);

% TC Dimensions
tc = tn; % Leg width
sl = 0.4*10^-3; % Spacing between unicouple legs
cs = 2*10^-3; % Spacing between Couples
cw = tc*2 +sl;
cvl = cw +cs;

F =4; % Couples Per Control Volume
E = floor(hxwC./(tc+(cs/10)));

if F>E
    fprintf('Number of Couples per Contol Volume too large\n')
    return
else
end

G = 2;
H  = floor(hxl./(cvl)); % Number of unicouples along fin direction

if G>H
    fprintf('Number of Couples per Contol Volume too large\n')
    return
else
end

C = G*F;

hxw = hxwC;
cvl = hxl; % Using TEBlock





% Input Checks

% PackFracColdSide = (tf_C*NC/hxwC)
% sC = (hxwC - (NC.*tf_C))./(NC-1)
if a ==1
    Tb = Tamb + Tcbg;
else
end

TempC = Tb+ Tcg; % [C]

[Q, P, V, I, Qmean ] = TEModule(Th,TempC, C);

QTE = Q;
% QTE = TEBlock(Th, TempC);

% Tb =40;

[Qc, hC, T] =HXMicrofinsCu(Tb, Tamb);

LegPackFrac = (ln*tn*2*G*F)/(hxlC*hxwC)
% Rth = (Tb-Tamb)/Qc;
  
%%

    while ( abs(Qc -QTE)> Err  )
        

         [Q, P, V, I, Qmean] = TEModule(Th, TempC, C);
         QTE = Q;
        [Qc, hC, TEMC] = HXMicrofinsCu(Tb, Tamb);%mean([Tcf(j), Tcf(j+1)]) );
       
        TempC = TEMC;
        if Qc -QTE> Err & Qc -QTE< (Err*DiffLim)   %Qc -QTE> Err & j ==1
          Tb = Tb -(TempCh);
        elseif Qc - QTE> (Err*DiffLim) 
            Tb = Tb -(ErrMult*TempCh);
        elseif (Qc -QTE) < (Err*-1) & Qc -QTE > (Err*- DiffLim)  %abs(Qc -QTE) > (Err*1) 
            Tb = Tb +(TempCh);
        elseif (Qc -QTE) < (Err*- DiffLim)
            Tb = Tb +(ErrMult*TempCh);
        end
        
%         
         
        diffcs  = Qc-QTE
        Tb
%         
%         progress = (j-1)/cv
    end
    
 progress = a/length(hnv)
    
 Rth = (TempC-Tamb)/mean([Qc,QTE]);
 Rth_TEM = (Th-TempC)/mean([Qc,QTE]);
 RthTEM = (Th -TempC)/Qmean;
 TotalP = P;
 Eff = (TotalP*100)/QTE
 
 
 Iv(a) =I;
 Vv(a) =V;
 Rthv(a) =Rth;
 Rth_TEMv(a) =Rth_TEM;
 RthTEMv(a) = RthTEM;
 TotalPv(a) = TotalP;
 Effv(a) =Eff;
 Qmeanv(a) = Qmean;
 
 
 end



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

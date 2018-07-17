% TEG Combined with Vertical Flat Plate Heat Sink
% Pasindu Gamarachchi - Email : pgamarachchi@gmail.com
clc
clear all

global hxl cvl NC tf_C Tcg Tcbg hxhC hxwC hxlC Cb  hn ln tn hp lp tp      

Th =200; % Hot Side Temperature
Tamb= 120;

hnv =  0.2*(10^-3): (0.1*(10^-3)): 2.0*(10^-3);
leng = length(hnv)
%%


%%Heat Sink Dimensions
NC = 7;
tf_C = 1.50*(10^-3);
hxhC = 0.15;  
hxwC = 0.04; 
hxl  = 0.04; 
hxlC = hxl;
Cb = 1*10^-2;
cvl = hxl; 
PackFracColdSide = (tf_C*NC/hxwC)
sC = (hxwC - (NC.*tf_C))./(NC-1)


% TC Dimensions
ln = 1.8*(10^-3);
tn = 1.8*(10^-3);
lp = 1.8*(10^-3);
tp = 1.8*(10^-3);
F =7;
G =7;
C = F*G;
FillFactor = (tn*ln*2*C)/(hxwC*hxlC)

% Guesses and Error
Err = 10^-6; 
Tcg = 1;
Tcbg =35; 
ErrMult = 3.5;
DiffLim = 800;  
ErrMult = 500;
TempCh = 0.5*Err;

if TempCh > Err
    fprintf('Unlikely to converge\n')
    return
end

% For each unicouple height
for i = 1:length(hnv)

    hn = hnv(i);
    hp = hn;

    if i ==1
        Tb = Tamb + Tcbg;
        TempC = Tb + Tcg;
    else
    end

    [QHT, QCT, P, n, Volt, I, Qm] = TEModuleDOE_Bi2Te3(Th,TempC,C);
    QTE = QCT;
    [Qc, Q_b, h_f, h_b ,TEMC ,TTip ,et_f, Qr] = HXFreeConvFlat(Tb, Tamb);

   
%% Iterative Condition
    while ( abs(Qc -QTE)> Err  )
        
        [QHT, QCT, P, n, Volt, I, Qm] = TEModuleDOE_Bi2Te3(Th,TempC,C);
        QTE = QCT;
        [Qc, Q_b, h_f, h_b ,TEMC ,TTip ,et_f, Qr] = HXFreeConvFlat(Tb, Tamb);
       
        TempC = TEMC;
        if Qc -QTE> Err & Qc -QTE< (Err*DiffLim)   
          Tb = Tb -(TempCh);
        elseif Qc - QTE> (Err*DiffLim) 
            Tb = Tb -(ErrMult*TempCh);
        elseif (Qc -QTE) < (Err*-1) & Qc -QTE > (Err*-DiffLim) 
            Tb = Tb +(TempCh);
        elseif (Qc -QTE) < (Err*-ErrMult)
            Tb = Tb +(ErrMult*TempCh);
        end
                 
        diffcs  = Qc-QTE
        Tb;
          
    end
 
    progress = i/length(hnv)

    TotalP = P;
    Eff = (TotalP*100)/QTE;
    TotalPv(i) = P;
    QHTv(i) = QHT;
    QCTv(i) = QCT;
    nv(i) =n;
    Voltv(i) = Volt;
    Iv(i) = I;
    Tbv(i) =Tb;
    TEMCv(i) = TEMC;

 
end



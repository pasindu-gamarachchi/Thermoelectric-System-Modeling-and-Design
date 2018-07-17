% Hot Side Heat Exchanger addeed to unicouple
% Separate function used for Hot Side Heat Exchanger and Unicouple
% Modified on : 07/18
clc
clear all
format long
tic


global hxl hxh hxw N Tin tf TEMW TEML tc sl cs cw E TempC cvl Thg F Tbg mf hn hp tn tp ln lp


% HX Dimensions & Input
hxl = 0.16; % [m] 
hxh = 0.02; % [m]
hxw = 0.04; % [m]

Nxv = [ 13, 14, 10,15               ]
tfv = [ 0.3*(10^-3), 0.3*(10^-3), 0.4*(10^-3), 0.4*(10^-3)]; %, 0.1*(10^-3), 0.1*(10^-3), 0.4*(10^-3), 0.4*(10^-3),  ]

for p = 1:length(Nxv)

N = Nxv(p); % Number of Fins

mf = 0.5*9.7*10^-3; % kg/s
Tin = 558; % [C]

% Fin Dimensions
% Material : Nickel
tf = tfv(p); % Fin thickness [m] 

% TEM Dimensions

TEMW = 27*10^-3; % [m]
TEML = 27*10^-3; % [m]

TempC = 94; % [C]

% Redudndant
tc = (1.5*10^-3); % Leg width
sl = 0.1*10^-3; % Spacing between unicouple legs
cs = 0.1*10^-3; % Spacing between Couples
cw = tc*2 +sl;
cvl = cw +cs;
% End of Redundancy

cvl = 2.5*(10^-3); %(40/6)*10^-3;
cv  = round(hxl./(cvl)); % Number of control volumes

F = 1; % Couples Per Control Volume

s = (hxw - (N.*tf))./(N+1); % Fin Spacing

E = 2*floor(TEMW./(cw + cs));

if F>E
    fprintf('Number of Couples per Contol Volume too large\n')
else
end

% Guesses and Error
Err = 0.005; % Convergence Error
Thg = 5; % Temperature Outlet Guess
Tbg = 180;  % Base/TEM Hot Side Temperature
TempCh = Err*0.9;
ErrMult = 5;
DiffLim = 1;

% Tvin = [ 400; 300; 200 ]
% for j=1:3
syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);
% 
knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

kCupf = int(5.25E-15*x^6 - 1.34E-11*x^5 + 1.26E-08*x^4 - 5.19E-06*x^3 + 8.71E-04*x^2 - 9.61E-02*x + 3.99E+02);
kCup = symfun(kCupf,x);


% Input Checks
PF = 100*(tf*N/hxw);
TotalCouples = F*cv;

%%
T = zeros(cv,1);
T(1) = Tin;

for j = 1:length(T)
     
   
    Tin = T(j);
    Th = Tin;
    Thi = Tin -Thg;
    Ti =Tin;

    cam = (ca(Th) - ca(Thi))./(Th-Thi);
    cam = double(cam);
  
    cp = cam;
  
    Tog = Ti - 10;

    if j==1;
        Tb = Ti-Tbg;
    else
    end
    
    Q = TEBlock(Tb,TempC);
%     [Q, P] = TEModule(Tb,TempC,1);
%     QT = Q*F;
    To = Ti - (Q/(mf*cp));
    Pf = 2*(cvl);
    Af = cvl*tf;
    
    Tfl = (Ti + To)/2;

%     knim = (kni(Tb) - kni(Tb-10))./(10);
%     knim = double(knim);
    
    kNim = (kni(Tb) -kni(Tb-10))./(10);
    kNim = double(kNim);

    % Fin Heat Transfer Calcs
    [h, v, Re, f_f, f_f2, prD, prD_O, prD_3] = convcoeffv3(T(j),mf);
    
%     h = 150;
     M = (sqrt(h*Pf*kNim*Af))*(Ti - Tb);
    
     Un_ar = (hxw- tf*N)*cvl;
     
    m = sqrt((h*Pf)/(kNim*Af));
    Qf = M*tanh(m*hxh/2)*((N-1)) + Un_ar*h*(Tb - Tfl) ;
%     [QTE,QCT, P, eff, Volt, I ,Qmean] = TEModule(Tb,TempC,1);
    QTET = TEBlock(Tb,TempC);

%     err = Tog- To;
    i=0;

    %%


    while abs(Qf-QTET)>Err
        
        Tb = Tb+(Err);
        [h, v, Re, f_f, f_f2, prD, prD_O, prD_3]= convcoeffv3(T(j),mf);
        M = (sqrt(h*Pf*kNim*Af))*((Ti- Tb));
        m = sqrt((h*Pf)/(kNim*Af));
        Qf = M*tanh(m*hxh/2)*((N-1)) + Un_ar*h*(Tb - Tfl) ;
        QTET = TEBlock(Tb,TempC);
        
        if Qf - QTET> Err & Qf -QTET< (Err*DiffLim)   %Qc -QTE> Err & j ==1
            Tb = Tb +(TempCh);
        elseif Qf - QTET> (Err*DiffLim) 
            Tb  = Tb  + (ErrMult*TempCh);
        elseif (Qf -QTET ) < (Err*-1) & Qf -QTET > (Err*- DiffLim)  %abs(Qc -QTE) > (Err*1) 
            Tb = Tb -(TempCh);
        elseif (Qf - QTET) < (Err*- DiffLim)
            Tb = Tb -(ErrMult*TempCh);
        end
        diff = Qf-QTET;


    end
    
    progress = j/length(T)
    QTETv(j) = QTET;
    
    hv(j) = h;
    vv(j) = v;
    Rev(j) = Re;
    prDv(j) = prD;
    f_fv(j) = f_f;
    f_f2v(j) = f_f2;
    prD_Ov(j) = prD_O;
    prDv(j) =prD;
    prD_3v(j) = prD_3;
%     Re2v(j) = Re2;
%     PTv(ja) =P*F;
    
%     cam = (ca(Th) - ca(Thi))./(Th-Thi);
%     cam = double(cam);
    
    To = Ti - (Qf./(mf*cp));
    Tbv(j) = Tb;
    Qhv(j) = Qf;
    T(j+1) = To;

   
end


%%
p = cv/4
Tbv = Tbv'
for i = 1:4
    
  
    TEMHMV(i) = mean( Tbv(((i-1)*p +1):p*i));
    QTEGHV(i) = sum( Qhv(((i-1)*p +1):p*i));
    
end

TEMHMV
Exh_DeltT = T(1) - T(end)
Avg_TEGQ = mean(QTEGHV)
Avg_TEMH = mean(TEMHMV)
% prD_tot = sum(prD_Ov);
prD_tot2 = sum(prDv);
prDO_tot = sum(prD_Ov);
prD3_tot = sum(prD_3v);
% fdum = fd*(10^6);
% fhmm = fh *(10^3);
tfmm = tf*(10^3);
dir = 'C:\Users\buddhimagamarachchi\Documents\Data\Matlab_Data\';
t = datetime;
d = datestr(t, 'dd-mmm-yyyy_HHMM');
name = ['tf_', num2str(tfmm), '_PF', num2str(PF), 'cv_ ', num2str(cv), '_ModelVerification_StrippedFinh_h'];
filename = [dir, name, d, '.mat'];
save(filename)


progress = p/length(Nxv)
end
% ln_cm = ln*100;
% tn_cm = tn*100;

%  progress = a/length(hnv)
% 
% end
 
% Pd = PTv./((hxl*100)*(hxw*100));
% RthHS = (Th - Tbv')./Qhv;
%  RthTEMv = (Tbv -TempC)./Qhv';
toc


% filename = [ 'PowerDensity_ThinFilm_', d_str, '.txt']
% save(filename, 'Pd', 'hnv')


% TC Dimensions
% a = (10*10^-6):(10*10^-6):(100*10^-6);
% b =(100*10^-6):(100*10^-6):(2000*10^-6);
% hnv  = [a,b];

% hnv = (10*10^-6):(10*10^-6): (20*10^-6)

% for a = 1:length(hnv)

% hn = hnv(a);
% hp = hn;
% tn = 1.2*(10^-3);
% ln = 1.2*(10^-3);
% tp = tn;
% lp = ln;


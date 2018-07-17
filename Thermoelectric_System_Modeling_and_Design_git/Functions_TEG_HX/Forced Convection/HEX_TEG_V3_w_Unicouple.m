% Hot Side Heat Exchanger addeed to unicouple
% Separate function used for Hot Side Heat Exchanger and Unicouple
% Modified on : 07/18
clc
clear all
format long
tic


global hxl hxh hxw N Tin tf TEMW TEML tc sl cs cw E TempC cvl Thg F Tbg mf hn hp tn tp ln lp
% HX Dimensions & Input
hxl = 2.406*(10^-3); % [m] 
hxh = 0.01; % [m]
hxw = 3.99*(10^-3); % [m]
N = 10; % Number of Fins

mf = 0.0005*9.7*10^-3; % kg/s
Tin = 350; % [C]

% Fin Dimensions

% Material : Nickel
tf = 0.10*10^-3; % Fin thickness [m] 

% TEM Dimensions

TEMW = 27*10^-3; % [m]
TEML = 27*10^-3; % [m]

 TempC = 40; % [C]

% TC Dimensions
a = (10*10^-6):(10*10^-6):(100*10^-6);
b =(100*10^-6):(100*10^-6):(2000*10^-6);
hnv  = [a,b];

% hnv = (10*10^-6):(10*10^-6): (20*10^-6)

for a = 1:length(hnv)

hn = hnv(a);
hp = hn;
tn = 1.2*(10^-3);
ln = 1.2*(10^-3);
tp = tn;
lp = ln;

tc = (1.5*10^-3); % Leg width
sl = 0.1*10^-3; % Spacing between unicouple legs
cs = 0.1*10^-3; % Spacing between Couples
cw = tc*2 +sl;
cvl = cw +cs;
% cvl = (40/6)*10^-3;
%cv  = floor(hxl./(cw+cs)); % Number of control volumes
cv  = round(hxl./(cvl)); % Number of control volumes

F = 1; % Couples Per Control Volume

s = (hxw - (N.*tf))./(N+1); % Fin Spacing
% Pf = 2*(cw+cs);
% Af = cv*tf;

E = 2*floor(TEMW./(cw + cs));

if F>E
    fprintf('Number of Couples per Contol Volume too large\n')
else
end

% Guesses and Error
Err = 0.001; % Convergence Error
Thg = 5; % Temperature Outlet Guess
Tbg = 258;  % Base/TEM Hot Side Temperature
TempCh = Err*0.9;
ErrMult = 5;
DiffLim = 1;

% Tvin = [ 400; 300; 200 ]
% for j=1:3
syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);
% 
% knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
% kni =  symfun(knif,x);

kCupf = int(5.25E-15*x^6 - 1.34E-11*x^5 + 1.26E-08*x^4 - 5.19E-06*x^3 + 8.71E-04*x^2 - 9.61E-02*x + 3.99E+02);
kCup = symfun(kCupf,x);


% Input Checks
PackFraction = (tf*N/hxw);
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

   

    % mf = 0.0097;
     cp = cam;

  
    Tog = Ti - 10;
    % Tbg = Ti- 1;
    if a==1;
        Tb = Ti-Tbg;
    else
    end
    
    [Q, P] = TEModule(Tb,TempC,1);
    QT = Q*F;
    To = Ti - (QT./(mf*cp));
    Pf = 2*(cvl);
    Af = cvl*tf;

%     knim = (kni(Tb) - kni(Tb-10))./(10);
%     knim = double(knim);
    
    kCum = (kCup(Tb) -kCup(Tb-10))./(10);
    kCum = double(kCum);

    % Fin Heat Transfer Calcs
    [heqn v Re] = convcoeff(T(j),mf);
    
    h = 150;
     M = (sqrt(h*Pf*kCum*Af))*(Ti - Tb);

    m = sqrt((h*Pf)/(kCum*Af));
    Qf = M*tanh(m*hxh/2)*((N-1));
    [QTE,QCT, P, eff, Volt, I ,Qmean] = TEModule(Tb,TempC,1);
    QTET = F*QTE;

%     err = Tog- To;
    i=0;

    %%


    while abs(Qf-QTET)>Err
        Tb = Tb+(Err);
        [heqn, v, Re]= convcoeff(T(j),mf);
        h = 150;
        M = (sqrt(h*Pf*kCum*Af))*((Ti- Tb));
        m = sqrt((h*Pf)/(kCum*Af));
        Qf = M*tanh(m*hxh/2)*((N-1));
        [QTE, QCT, P, eff,Volt,I, Qmean] = TEModule(Tb,TempC,1);
        QTET =QTE*F;
        
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
%         if Qf<QTET
%             fprintf('Tbg too small\n');
%         else
%         end
%         i=i+1;

    end
    QTETv(j,a) = QTET;
    
    hv(j,a) = h;
    vv(j,a) = v;
    Rev(j,a) = Re;
    PTv(j,a) =P*F;
    
%     cam = (ca(Th) - ca(Thi))./(Th-Thi);
%     cam = double(cam);
    
    To = Ti - (Qf./(mf*cp));
    Tbv(j,a) = Tb;
    Qhv(j,a) = Qf;
    T(j+1,a) = To;

   
end

% ln_cm = ln*100;
% tn_cm = tn*100;

 progress = a/length(hnv)

end
 
Pd = PTv./((hxl*100)*(hxw*100));
% RthHS = (Th - Tbv')./Qhv;
%  RthTEMv = (Tbv -TempC)./Qhv';
toc

d_str = date;
filename = [ 'PowerDensity_ThinFilm_', d_str, '.txt']
save(filename, 'Pd', 'hnv')



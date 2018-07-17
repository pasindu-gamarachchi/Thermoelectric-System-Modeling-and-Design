% Hot Side Heat Exchanger and TE Unicouple model combined
% Pasindu Gamarachchi - Email: pgamarachchi@gmail.com
clc
clear all


global hxl hxh hxw N Tin tf TEMW TEML TempC cvl Thg  Tbg mf 

% HX Dimensions & Input
hxl = 0.16; 
hxh = 0.02; 
hxw = 0.04; 
N = 30; 
mf = 0.5*9.7*10^-3; 
Tin = 558; 
tf = 0.1*10^-3; 
s = (hxw - (N.*tf))./(N+1); % Fin Spacing

% TEM Dimensions

TEMW = 40*10^-3; % [m]
TEML = 40*10^-3; % [m]

TempC = 94; % [C]

cvl = 40*(10^-3); 
cv  = round(hxl./(cvl)); % Number of control volumes

% Guesses and Error
Err = 0.005; % Convergence Error
Thg = 5; % Temperature Outlet Guess
Tbg = 180;  % Base/TEM Hot Side Temperature
Tdg = 10;
TempCh = Err*0.9;
ErrMult = 5;
DiffLim = 1;

% Fin Material and Fluid Properties
syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);
% 
knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

%%
T = zeros(cv,1);
T(1) = Tin;

% Calculations for each control volume
for j = 1:length(T)
     
   
    Tin = T(j);
    Th = Tin;
    Thi = Tin -Thg;
    Ti =Tin;

    cam = (ca(Th) - ca(Thi))./(Th-Thi);
    cam = double(cam);
  
    cp = cam;
  
    Tog = Ti - Tdg;

    if j==1;
        Tb = Ti-Tbg;
    else
    end
    
    [Q,P,Eff] = TEModule(Tb,TempC, Un);% Function with finite element unicouple function and number of unicouples

    To = Ti - (Q/(mf*cp));
    Pf = 2*(cvl);
    Af = cvl*tf;
    
    Tfl = (Ti + To)/2;
    
    kNim = (kni(Tb) -kni(Tb-10))./(10);
    kNim = double(kNim);

    % Fin Heat Transfer Calcs
    [h, v, Re] = convcoeff(Tfl,mf); % Function with Duct Convection Coefficient
    
    M = (sqrt(h*Pf*kNim*Af))*(Ti - Tb);
    Un_ar = (hxw- tf*N)*cvl;
     
    m = sqrt((h*Pf)/(kNim*Af));
    Qf = M*tanh(m*hxh/2)*((N-1)) + Un_ar*h*(Tb - Tfl) ;

    [QTET, P, Eff] = TEModule(Tb,TempC, Un);

    i=0;
    
    % Convergence Requirement
    while abs(Qf-QTET)>Err
        
        Tb = Tb+(Err);
        [h, v, Re]= convcoeffv3(T(j),mf);
        M = (sqrt(h*Pf*kNim*Af))*((Ti- Tb));
        m = sqrt((h*Pf)/(kNim*Af));
        Qf = M*tanh(m*hxh/2)*((N-1)) + Un_ar*h*(Tb - Tfl) ;
        [QTET, P, Eff] = TEModule(Tb,TempC, Un);
        
        if Qf - QTET> Err & Qf -QTET< (Err*DiffLim)  
            Tb = Tb +(TempCh);
        elseif Qf - QTET> (Err*DiffLim) 
            Tb  = Tb  + (ErrMult*TempCh);
        elseif (Qf -QTET ) < (Err*-1) & Qf -QTET > (Err*- DiffLim) 
            Tb = Tb -(TempCh);
        elseif (Qf - QTET) < (Err*- DiffLim)
            Tb = Tb -(ErrMult*TempCh);
        end
        diff = Qf-QTET;

    end
    
    progress = j/length(T)
    QTETv(j) = QTET;
    Pv (j) = P;
    Effv(j) = Eff;
    
    hv(j) = h;
    vv(j) = v;
    Rev(j) = Re;
  
    To = Ti - (Qf./(mf*cp));
    Tbv(j) = Tb;
    Qhv(j) = Qf;
    T(j+1) = To;
   
end





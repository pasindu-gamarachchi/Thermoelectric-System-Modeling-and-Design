% 
clc
clear all
format long
tic


global hxl hxh hxw N Tin tf TEMW TEML tc sl cs cw E TempC cvl Thg F
% HX Dimensions & Input
hxl = 3.75*(10^-3); % [m]
hxh = 0.02; % [m]
hxw = 0.04; % [m]
N = 40; % Number of Fins

mf = 0.5*9.7*10^-3; % kg/s
Tin = 558; % [C]

% Fin Dimensions

% Material : Nickel
tf = 0.1*10^-3; % Fin thickness [m] 

% TEM Dimensions

TEMW = 27*10^-3; % [m]
TEML = 27*10^-3; % [m]

 TempC = 94; % [C]

% TC Dimensions
tc = (1.8*10^-3); % Leg width
sl = 0.4*10^-3; % Spacing between unicouple legs
cs = 3*10^-3; % Spacing between Couples
cw = tc*2 +sl;
cvl = cw +cs;
% cvl = (40/6)*10^-3;
%cv  = floor(hxl./(cw+cs)); % Number of control volumes
cv  = floor(hxl./(cvl)); % Number of control volumes

F = 4; % Couples Per Control Volume

s = (hxw - (N.*tf))./(N+1); % Fin Spacing
% Pf = 2*(cw+cs);
% Af = cv*tf;

E = 2*floor(TEMW./(cw + cs));

if F>E
    fprintf('Number of Couples per Contol Volume too large\n')
else
end

% Guesses and Error
Err = 0.05; % Convergence Error
Thg = 5; % Temperature Outlet Guess
Tbg = 48;  % Base/TEM Hot Side Temperature

% Tvin = [ 400; 300; 200 ]
% for j=1:3
syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

% Input Checks
PackFraction = (tf*N/hxw)
TotalCouples = F*cv

%%
T = zeros(cv,1);
T(1) = Tin;

for j = 1:length(T)
     
   
    Tin = T(j)
    Th = Tin;
    Thi = Tin -Thg;
    Ti =Tin;




    cam = (ca(Th) - ca(Thi))./(Th-Thi);
    cam = double(cam);

   

    % mf = 0.0097;
     cp = cam;

  
    Tog = Ti - 10;
    % Tbg = Ti- 1;
    Tb = Ti-Tbg;
    
    [Q P] = Unicouple(Tb,TempC);
    QT = Q*F;
    To = Ti - (QT./(mf*cp));
    Pf = 2*(cvl);
    Af = cvl*tf;

    knim = (kni(Tb) - kni(Tb-10))./(10);
    knim = double(knim);

    % Fin Heat Transfer Calcs
    [h v Re] = convcoeff(T(j),mf);
     M = (sqrt(h*Pf*knim*Af))*(Tb - Ti);

    m = sqrt((h*Pf)/(knim*Af));
    Qf = M*tanh(m*hxh/2)*((N-1));
    [QTE P] = Unicouple(Tb,TempC);
    QTET = F*QTE;

%     err = Tog- To;
    i=0;

    %%


    while abs(Qf-QTET)>Err
        Tb = Tb+(Err);
        [h v Re]= convcoeff(T(j),mf);
        M = (sqrt(h*Pf*knim*Af))*((mean([Ti,To])- Tb));
        m = sqrt((h*Pf)/(knim*Af));
        Qf = M*tanh(m*hxh/2)*((N-1));
        [QTE P] = Unicouple(Tb,TempC);
        QTET =QTE*F;
        if Qf<QTET
            fprintf('Tbg too small\n');
        else
        end
        i=i+1;

    end
    QTETv(j) = QTET;
    
    hv(j) = h;
    vv(j) = v;
    Rev(j) = Re;
    PTv(j) =P*F;
    
%     cam = (ca(Th) - ca(Thi))./(Th-Thi);
%     cam = double(cam);
    
    To = Ti - (Qf./(mf*cp));
    Tbv(j) = Tb;
    Qhv(j) = Qf;
    T(j+1) = To;

    progress = j/cv
end


% end
toc
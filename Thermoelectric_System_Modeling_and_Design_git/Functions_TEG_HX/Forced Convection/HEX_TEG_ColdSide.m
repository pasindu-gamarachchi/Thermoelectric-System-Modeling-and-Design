% 
clc
clear all
format long
tic

% New Version 03/23 
% Cold Side Heat Exchanger Added
% Separate Function Used for Cold Side Heat Exchanger




global hxl hxh hxw N Tin tf TEMW TEML tc sl cs cw E TempC cvl Thg Tcg mfC Tcbg NC tf_C
% HX Dimensions & Input
hxl = 0.16; % [m]
hxh = 0.02; % [m]
hxw = 0.04; % [m]
N = 30; % Number of Fins

mf = 0.5*9.7*10^-3; % kg/s
Tin = 558; % [C]

% Fin Dimensions

% Material : Nickel
tf = 0.4*10^-3; % Fin thickness [m] 

% TEM Dimensions

TEMW = 27*10^-3; % [m]
TEML = 27*10^-3; % [m]

 TempC = 48; % [C]

% TC Dimensions
tc = (1.8*10^-3); % Leg width
sl = 0.4*10^-3; % Spacing between unicouple legs
cs = 2*10^-3; % Spacing between Couples
cw = tc*2 +sl;
%cvl = cw +cs;
cvl = (5)*10^-3; % Using TEBlock
%cv  = floor(hxl./(cw+cs)); % Number of control volumes
cv  = floor(hxl./(cvl)); % Number of control volumes



s = (hxw - (N.*tf))./(N+1); % Fin Spacing
% Pf = 2*(cw+cs);
% Af = cv*tf;

E = 2*floor(TEMW./(cw + cs));
%

% Cold Side Input
TCoIn= 20;
Tcg  = 8;
mfC  = 0.5; % kg/s

Tcbg= 70;
NC = 20;
tf_C= 0.1*10^-3; 
cvl = 5*10^-3;


% Guesses and Error
Err = 0.1; % Convergence Error
Thg = 5; % Temperature Outlet Guess
Tbg = 48;  % Base/TEM Hot Side Temperature
Tcfg = 5;


% Tvin = [ 400; 300; 200 ]
% for j=1:3
syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);

% Input Checks
PackFraction_Cold = (tf_C*NC/hxw);

%%
T = zeros(cv,1);
T(1) = Tin;
Tcf = zeros(cv,1);
Tcf(1) = TCoIn;
TEMC = zeros(cv,1);
Qcv = zeros(cv,1);
% Tb = Tin-Tbg;
for j = 1:length(T)
     
    Tb =500
%     Tin = T(j)
%     Th = Tin;
%     Thi = Tin -Thg;
%     Ti =Tin;
    Tcf(j+1) = Tcf(j) + Tcfg
%     Tb = Tb -5;


%     cam = (ca(Th) - ca(Thi))./(Th-Thi);
%     cam = double(cam);

   

    % mf = 0.0097;
%      cp = cam;

% 
%     Tog = Ti - 10;
    % Tbg = Ti- 1;
%     Tb = Ti-Tbg;

     Q = TEBlock(Tb,TempC)
%     To = Ti - (Q./(mf*cp));
    Pf = 2*(cvl)
    Af = cvl*tf

    knim = (kni(Tb) - kni(Tb-10))./(10);
    knim = double(knim);

    % Fin Heat Transfer Calcs
    %[h v Re] = convcoeff(T(j),mf);
     %M = (sqrt(h*Pf*knim*Af))*(Tb - Ti);

    %m = sqrt((h*Pf)/(knim*Af));
    %Qf = M*tanh(m*hxh/2)*((N-1));
    QTE = TEBlock(Tb, Tcf(j));


%     err = Tog- To;
    i=0;
    Qc = HXCold(TempC, Tcf(j));
    %%
   

    while abs(Qc -QTE)> Err 
        
%         TempC = TempC -0.1;
        %[h v Re]= convcoeff(T(j),mf);
        %M = (sqrt(h*Pf*knim*Af))*((Ti- Tb));
        %m = sqrt((h*Pf)/(knim*Af));
        %Qf = M*tanh(m*hxh/2)*((N-1));
        QTE = TEBlock(500, TempC)
        [Qc ToutC h] = HXCold(TempC, mean([Tcf(j), Tcf(j+1)])  );
       
        if Qc -QTE> Err
          TempC = TempC -(Err*0.2)
        elseif Qc -QTE <-Err
            TempC = TempC +(Err*0.2)
        
        end
        
%         if (Qf-QTE)>0.1
%           Tb = Tb+0.02;
%         elseif (Qf- QTE) <0.1
%             Tb = Tb -0.02;
% %         end
%         
%                 
%         
%         if ( abs(Qc -QTE) +abs(Qf-QTE) ) < 0.21
%             break
%         else
%         end
%         if Qc<0
%             break
%             fprintf('Error\n')
%         else
%         end
%          diffhs = Qf-QTE
        diffcs  = Qc-QTE
%         progress = (j-1)/cv
    end
    QTEv(j) = QTE;
    
%     hv(j) = h;
%     vv(j) = v;
%     Rev(j) = Re;
    TEMC(j) = TempC
    Tcf(j+1) = ToutC
    hcv(j) =h;
   % To = Ti - (Qf./(mf*cp));
%     Tbv(j) = Tb;
%     Qhv(j) = Qf;
    Qcv(j) = Qc;
%     T(j+1) = To;
    

    progress = j/cv
end
% end
toc
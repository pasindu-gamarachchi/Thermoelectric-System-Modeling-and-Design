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


global hxl  TempC cvl  NC tf_C  Tcg Tcbg hxhC hxwC hxlC Cb hxw


Tamb= 20;


% Cold Side Dimensions
NC = 7;
tf_C = 0.214*(10^-3); 

hxhC = 0.02;  %m
hxwC = 0.04; %m
hxl  = 0.04; %m
hxlC = hxl;
Cb = 0.5*10^-2; % m Cold Base Thickness

% TC Dimensions
tc = (1.8*10^-3); % Leg width
sl = 0.4*10^-3; % Spacing between unicouple legs
cs = 1*10^-3; % Spacing between Couples
cw = tc*2 +sl;
cvl = cw +cs;

hxw = hxwC;
cvl = hxl; % Using TEBlock
cv  = floor(hxl./(cvl)); % Number of control volumes

Th = 37; % Hot Side Temperature

% Guesses and Error
Err = 0.015; % Convergence Error

Tcg = 0.004; % [C] Guess for temperature of TEMC
Tcbg = 12; %  [C] Guess for temperature of Fin Base
ErrMult = 5.5;



syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);



% Input Checks

PackFracColdSide = (tf_C*NC/hxwC)
sC = (hxwC - (NC.*tf_C))./(NC-1)
Tb = Tamb + Tcbg;
TempC = Tb+0.04; % [C]

QTE = TEBlock(Th, TempC);

% Tb =40;

[Qc hC] = HXFreeConv(Tb, Tamb)

Rth = (Tb-Tamb)/Qc
    %%


    while ( abs(Qc -QTE)> Err  )
        

         QTE = TEBlock(Th, TempC);
        [Qc hC TEMC Ra] = HXFreeConv(Tb, Tamb);%mean([Tcf(j), Tcf(j+1)]) );
       
        TempC = TEMC;
        if Qc -QTE> Err & Qc -QTE< (Err*ErrMult)   %Qc -QTE> Err & j ==1
          Tb = Tb -(Err*0.2);
        elseif Qc - QTE> (Err*ErrMult) 
            Tb = Tb -(Err*8.0);
        elseif (Qc -QTE) < (Err*-1) & Qc -QTE > (Err*-ErrMult)  %abs(Qc -QTE) > (Err*1) 
            Tb = Tb +(Err*0.2);
        elseif (Qc -QTE) < (Err*-ErrMult)
            Tb = Tb +(Err*8.0);
        end
        
%         
         
        diffcs  = Qc-QTE
        Tb;
%         
%         progress = (j-1)/cv
    end
 
 %%
 Rth = (Tb-Tamb)/Qc;
 
%  fprintf('Thermal Resistance is %
 fprintf('Thermal Resistance is %d C/W\n', Rth)
fprintf('TEM Cold Side Temp is %d C\n', TEMC)
fprintf('Cold Side Fin Base Temp is %d C\n', Tb)
fprintf('Heat Transferred is %d W\n', Qc)
fprintf('Convection Coeff is %d W/m^2-C\n', hC)
fprintf('Fin Spacing is %d mm\n', sC)

    
% hxh hxw N Tin tf TEMW TEML tc sl cs cw E

% HX Dimensions & Input
% hxl = 0.16; % [m]
% hxh = 0.02; % [m]
% hxw = 0.04; % [m]
% N = 30; % Number of Fins

% mf = 0.5*9.7*10^-3; % kg/s
% Tin = 558; % [C]

% Fin Dimensions

% Material : Nickel
% tf = 0.2*10^-3; % Fin thickness [m] 

% TEM Dimensions

% TEMW = 27*10^-3; % [m]
% TEML = 27*10^-3; % [m]

% Cold Side Input

% mfC  =0.04; % kg/s

% TC Dimensions
% tc = (1.8*10^-3); % Leg width
% sl = 0.4*10^-3; % Spacing between unicouple legs
% cs = 2*10^-3; % Spacing between Couples
% cw = tc*2 +sl;
%cvl = cw +cs;

%cv  = floor(hxl./(cw+cs)); % Number of control volumes

% s = (hxw - (N.*tf))./(N+1); % Fin Spacing
% Pf = 2*(cw+cs);
% Af = cv*tf;

% E = 2*floor(TEMW./(cw + cs));
%

% Tvin = [ 400; 300; 200 ]
% for j=1:3
% PackFracHotSide = (tf*N/hxw);

% Initialize Temperature Vectors
% T = zeros(cv,1);
% Tcf = zeros(cv,1);
% TEMC = zeros(cv,1);
% Qcv = zeros(cv,1);

%%

% T(1) = Tin;
% Tcf(1) = TCoIn;
% Tb = Tin-Tbg;

% Loop Independent Calculations
% Pf = 2*(cvl);
% Af = cvl*tf;

% for j = 1:length(T)
%      
%    
%     Ti = T(j);
% %     Th = Tin;
%     Thi = Ti -Thg;
% %     Ti =Ti;
%     
%     if j~=1
% %          Tb = Tb -5;
% %         TempC = TempC +1;
%     else
%     end
    

    

%     cam = (ca(Ti) - ca(Thi))./(Ti-Thi);
%     cam = double(cam);
% 
%      cp = cam;


%     Tog = Ti - 10;


%     Q = TEBlock(Th,TempC);
%     To = Ti - (Q./(mf*cp));


%     knim = (kni(Tb) - kni(Tb-10))./(10);
%     knim = double(knim);

    % Fin Heat Transfer Calcs
%     [h v Re] = convcoeff(T(j),mf);
%      M = (sqrt(h*Pf*knim*Af))*(Tb - Ti);
% 
%     m = sqrt((h*Pf)/(knim*Af));
%     Qf = M*tanh(m*hxh/2)*((N-1));

%     err = Tog- To;
%     i=0;

% Tcg = 0.5; % Cold Side fluid delta T guess

%         TempC = TempC -0.1;
%         [h v Re]= convcoeff(T(j),mf);
%         M = (sqrt(h*Pf*knim*Af))*((Ti- Tb));
%         m = sqrt((h*Pf)/(knim*Af));
%         Qf = M*tanh(m*hxh/2)*((N-1));


% if (Qf-QTE)>Err & Qf -QTE< (Err*ErrMult)   %& j ==1   % (Qf-QTE)>Err  & j ==1
%           Tb = Tb+(Err*0.2);
%         elseif Qf - QTE> (Err*ErrMult) 
%             Tb = Tb +(Err*10.0);
%         elseif  (Qf- QTE)< -Err & (Qf-QTE) >(Err*-ErrMult)
%             Tb = Tb -(Err*0.1);
%         elseif (Qf -QTE) <(-Err*ErrMult)
%             Tb  = Tb -(Err*10.0);
%         end
        
                
%         
%         if ( abs(Qc -QTE) +abs(Qf-QTE) ) < Err*2.0
%             break
%         else
%         end
%         if Qc<0
%             break
%             fprintf('Error\n')
%         else
%         end
        
%         Qf
%         QTE
%         Qc
%          diffhs = Qf-QTE
%         Tb;
% Thg = 45; % Temperature Outlet Guess
% Tbg =3;  % Base/TEM Hot Side Temperature
% Tcfg= 0.5; 
%     hCv(j) = hC;
%     vCv(j) = vC;
%     QTEv(j) = QTE;
%     
%     hv(j) = h;
%     vv(j) = v;
%     Rev(j) = Re;
%     TEMC(j) = TempC;
%     Tcf(j+1) = ToutC;
%     
%     if j ==1
%         vC
%         hC
%     else
%     end

%     To = Ti - (Qf./(mf*cp));
%     Tbv(j) = Tb;
%     Qhv(j) = Qf;
%     Qcv(j) = Qc;
%     T(j+1) = To;
%     Recv(j)= ReC;
%     
% 
%     progress = j/cv
% end
% % end
% toc


%%
% x =linspace(0,hxl,cv);
% y = linspace(0,hxl, cv+1);
% 
% figure (1)
% plot(x,TEMC, 'b--')
% hold on
% plot(x, Tbv, 'r--')
% plot(y,T, 'g--')
% plot(y,Tcf, 'k--')
% xlabel('Position along x [m]')
% ylabel('Temperature [C]')
% 
% figure (2)
% semilogy(x,hCv, 'b--')
% hold on
% semilogy(x, hv, 'r--')
% xlabel('Position along x [m]')
% ylabel('Convection Coefficient [W/m^2-K]')
% 
% HotFluid =mean(T)
% TEMHot = mean(Tbv)
% TEMCold = mean(TEMC)
% ColdFluid = mean(Tcf)
% ColdConvection = mean(hCv)
% ColdReynolds = mean(Recv)
% ColdVelocity= mean(vCv)
% ColdSideBiot = ColdConvection*tf_C/(235)
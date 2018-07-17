% 
clc
clear 
close all
format long
tic

% New Version 03/23 
% Cold Side Heat Exchanger Added
% Separate Function Used for Cold Side Heat Exchanger
% Changes made 04/21, Better programming protocals implemented
% Changes made 04/28, unneccessary variabls removed



global hxl hxh hxw N Tin tf TEMW TEML tc sl cs cw E TempC cvl Thg NC tf_C mfC Tcg Tcbg hxhC
% HX Dimensions & Input
hxl = 0.16; % [m]
hxh = 0.02; % [m]
hxw = 0.04; % [m]
N = 30; % Number of Fins

mf = 0.5*9.7*10^-3; % kg/s
Tin = 558; % [C]

% Fin Dimensions

% Material : Nickel
tf = 0.2*10^-3; % Fin thickness [m] 

% TEM Dimensions

TEMW = 27*10^-3; % [m]
TEML = 27*10^-3; % [m]



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
TCoIn= 82;

NC = 8;
tf_C= 1*(10^-3); 
mfC  =0.04; % kg/s


hxhC = 1*10^-3; 

% Guesses and Error
Err = 1; % Convergence Error
Thg = 45; % Temperature Outlet Guess
Tbg =3;  % Base/TEM Hot Side Temperature
% Tcfg= 0.5; 
Tcg = 0.5; % Cold Side fluid delta T guess
Tcbg = 10;
ErrMult = 5.5;


% Tvin = [ 400; 300; 200 ]
% for j=1:3
syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01); % Thermal Conductivity of Nickel
kni =  symfun(knif,x);



% Input Checks
PackFracHotSide = (tf*N/hxw);
PackFracColdSide = (tf_C*NC/hxw);
 TempC = TCoIn + Tcbg; % [C]
% Initialize Temperature Vectors
T = zeros(cv,1);
Tcf = zeros(cv,1);
TEMC = zeros(cv,1);
Qcv = zeros(cv,1);

%%

T(1) = Tin;
Tcf(1) = TCoIn;
Tb = Tin-Tbg;

% Loop Independent Calculations
Pf = 2*(cvl);
Af = cvl*tf;

for j = 1:length(T)
     
   
    Ti = T(j);
%     Th = Tin;
    Thi = Ti -Thg;
%     Ti =Ti;
    
    if j~=1
%          Tb = Tb -5;
%         TempC = TempC +1;
    else
    end
    

    

    cam = (ca(Ti) - ca(Thi))./(Ti-Thi);
    cam = double(cam);

     cp = cam;


    Tog = Ti - 10;


    Q = TEBlock(Tb,TempC);
    To = Ti - (Q./(mf*cp));


    knim = (kni(Tb) - kni(Tb-10))./(10);
    knim = double(knim);

    % Fin Heat Transfer Calcs
    [h v Re] = convcoeff(T(j),mf);
     M = (sqrt(h*Pf*knim*Af))*(Tb - Ti);

    m = sqrt((h*Pf)/(knim*Af));
    Qf = M*tanh(m*hxh/2)*((N-1));
    QTE = TEBlock(Tb, Tcf(j));


%     err = Tog- To;
    i=0;
    Qc = HXColdWater(TempC, Tcf(j));
    %%


    while ( abs(Qc -Qf)> Err || abs(Qf-QTE)> Err  )
        
%         TempC = TempC -0.1;
        [h v Re]= convcoeff(T(j),mf);
        M = (sqrt(h*Pf*knim*Af))*((Ti- Tb));
        m = sqrt((h*Pf)/(knim*Af));
        Qf = M*tanh(m*hxh/2)*((N-1));
        QTE = TEBlock(Tb, TempC);
        [Qc ToutC hC vC ReC] = HXColdWater(TempC, Tcf(j));%mean([Tcf(j), Tcf(j+1)]) );
       
        if Qc -QTE> Err & Qc -QTE< (Err*ErrMult)   %Qc -QTE> Err & j ==1
          TempC = TempC -(Err*0.2);
        elseif Qc - QTE> (Err*ErrMult) 
            TempC = TempC -(Err*4.0);
        elseif (Qc -QTE) < (Err*-1) & Qc -QTE > (Err*-ErrMult)  %abs(Qc -QTE) > (Err*1) 
            TempC = TempC +(Err*0.2);
        elseif (Qc -QTE) < (Err*-ErrMult)
            TempC = TempC +(Err*4.0);
        end
        
        if (Qf-QTE)>Err & Qf -QTE< (Err*ErrMult)   %& j ==1   % (Qf-QTE)>Err  & j ==1
          Tb = Tb+(Err*0.2);
        elseif Qf - QTE> (Err*ErrMult) 
            Tb = Tb +(Err*10.0);
        elseif  (Qf- QTE)< -Err & (Qf-QTE) >(Err*-ErrMult)
            Tb = Tb -(Err*0.1);
        elseif (Qf -QTE) <(-Err*ErrMult)
            Tb  = Tb -(Err*10.0);
        end
        
                
        
        if ( abs(Qc -QTE) +abs(Qf-QTE) ) < Err*2.0
            break
        else
        end
        if Qc<0
            break
            fprintf('Error\n')
        else
        end
        
%         Qf
%         QTE
%         Qc
         diffhs = Qf-QTE
        Tb;
         
        diffcs  = Qc-QTE
        TempC;
%         
%         progress = (j-1)/cv
    end
    
    hCv(j) = hC;
    vCv(j) = vC;
    QTEv(j) = QTE;
    
    hv(j) = h;
    vv(j) = v;
    Rev(j) = Re;
    TEMC(j) = TempC;
    Tcf(j+1) = ToutC;
    
    if j ==1
        vC
        hC
    else
    end

    To = Ti - (Qf./(mf*cp));
    Tbv(j) = Tb;
    Qhv(j) = Qf;
    Qcv(j) = Qc;
    T(j+1) = To;
    Recv(j)= ReC;
    

    progress = j/cv
end
% end
toc

%%
x =linspace(0,hxl,cv);
y = linspace(0,hxl, cv+1);

figure (1)
plot(x,TEMC, 'b--')
hold on
plot(x, Tbv, 'r--')
plot(y,T, 'g--')
plot(y,Tcf, 'k--')
xlabel('Position along x [m]')
ylabel('Temperature [C]')

figure (2)
semilogy(x,hCv, 'b--')
hold on
semilogy(x, hv, 'r--')
xlabel('Position along x [m]')
ylabel('Convection Coefficient [W/m^2-K]')

HotFluid =mean(T)
TEMHot = mean(Tbv)
TEMCold = mean(TEMC)
ColdFluid = mean(Tcf)
ColdConvection = mean(hCv)
ColdReynolds = mean(Recv)
ColdVelocity= mean(vCv)
ColdSideBiot = ColdConvection*tf_C/(235)
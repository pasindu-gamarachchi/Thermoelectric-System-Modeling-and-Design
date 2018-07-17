%Iterator
% clc
% clear all


% function [To Q] =Iterator(Ti,mf)


% global hxw hxl tc sl cs cw E TempC



% TEMW = 27*10^-3; % [m]
% TEML = 27*10^-3; % [m]
% 
%  Ti = 558;
%  hxw = 0.04;
%  tc = 1.8*10^-3;
%  TempC = 94;
% N = 5;
% h = 188;
% P = 0.0042;
% A = 9*10^-5;
% k = 20;
% % sl = 0.4*10^-3; % Spacing between unicouple legs
% L = 3.6*10^-3;
% tc = 1.8*10^-3; 
% O = L/(tc.*2);% Couples per control volume along HX
% cw = tc*2 +sl;
% cs = 2*10^-3;
% E = 2*round(hxw./(cw+cs)); % Total number of couples for the cv, 2 to account for top and bottom
Th = Tin;
Thi = Tin -10;
Ti =Tin;

syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01);
kni =  symfun(knif,x);

 




cam = (ca(Th) - ca(Thi))./(Th-Thi);
cam = double(cam);

knim = (kni(Th) - kni(Thi))./(Th-Thi);
knim = double(knim);

% mf = 0.0097;
 cp = cam;

% TC Dimensions
% tc = (1.8*10^-3); % Leg width
% sl = 0.4*10^-3; % Spacing between unicouple legs
% cs = 2*10^-3; % Spacing between Couples
% cw = tc*2 +sl;
% cv  = floor(hxl./(cw+cs)); % Number of control volumes

% E = 2*floor(TEMW./(cw + cs));
Tog = Ti - 40;
% Tbg = Ti- 1;
Tb = Ti-50;
% m = sqrt((h.*P)./(k.*A));
% M = m.*((Ti - Tb));
% Q = 2.*N.*O.*M.*tanh(m.*L);
% Q = Unicouple(Tb,TempC);
 Q = TEBlock(Tb,TempC)
To = Ti - (Q./(mf*cp));



% Fin Heat Transfer Calcs
h = convcoeff(Th,mf);
M = sqrt(h*Pf*knim*Af)*(Tb - Ti);
m = sqrt((h*Pf)/(knim*Af));
Qf = -M*tanh(m*hxh/2)*(2*(N-1))
QTE = TEBlock(Tb, TempC)


err = Tog- To;
i=0;

%%


while abs(Qf-QTE)>1
    Tb = Ti 


%     
%     while abs(Tog -To) > 1
% 
%     Tog = Tog +1;
%     Tb  = (Ti+Tog)./2;
%     Q = TEBlock(Tb,TempC); % E*Unicouple(Tb,94);
%     % m = sqrt((h.*P)./(k.*A));
%     % M = m.*(Ti - Tb);
%     % Q = 2.*N.*O*M.*tanh(m.*L);
%     To = Ti - (Q./(mf*cp));
%     % Tog = Tog +0.1;
%     i= i+1;
%     end
%  i;
%  To;

% Qu = Unicouple(Tbg, 100)
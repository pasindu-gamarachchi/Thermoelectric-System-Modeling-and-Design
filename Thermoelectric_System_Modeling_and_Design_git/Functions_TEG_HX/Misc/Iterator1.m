%Iterator
% clc
% clear all

%% This function is used to estimate the initial guess
% function [To Qf] =Iterator1(Tin,mf)
% 
% 
% global hxw hxl tc sl cs cw E TempC cv hxh



Th = Tin;
Thi = Tin -5;
Ti =Tin;

syms x
caf = int( 3.134242E-10*x^4 - 8.519344E-07*x^3 + 7.480582E-04*x^2 - 3.006360E-02*x + 1.007301E+03); % specific Heat of air
ca =  symfun(caf,x);

knif = int( -9.32400932E-11*x^4 + 1.13247863E-07*x^3 + 6.33449883E-05*x^2 -9.47163947E-02*x + 8.13811189E+01);
kni =  symfun(knif,x);

 




cam = (ca(Th) - ca(Thi))./(Th-Thi);
cam = double(cam);



% mf = 0.0097;
 cp = cam;

% TC Dimensions
% tc = (1.8*10^-3); % Leg width
% sl = 0.4*10^-3; % Spacing between unicouple legs
% cs = 2*10^-3; % Spacing between Couples
% cw = tc*2 +sl;
% cv  = floor(hxl./(cw+cs)); % Number of control volumes

% E = 2*floor(TEMW./(cw + cs));
Tog = Ti - 10;
% Tbg = Ti- 1;
Tb = Ti-65;

knim = (kni(Tb) - kni(Tb-10))./(10);
knim = double(knim);
% m = sqrt((h.*P)./(k.*A));
% M = m.*((Ti - Tb));
% Q = 2.*N.*O.*M.*tanh(m.*L);
[Q P]= Unicouple(Tb,TempC);
%  Q = TEBlock(Tb,TempC);
To = Ti - (4*Q./(mf*cp));
Pf = 2*(cvl);
Af = cvl*tf;


% Fin Heat Transfer Calcs
h = convcoeff(Th,mf);
 M = (sqrt(h*Pf*knim*Af))*(Tb - Ti);

m = sqrt((h*Pf)/(knim*Af));
Qf = -M*tanh(m*hxh/2)*((N-1));
[QTE P]= Unicouple(Tb,TempC);

QTET =F*QTE;
err = Tog- To;
i=0;

%%


while abs(Qf-QTET)>1
    Tb = Tb+1;
    [h v] = convcoeff(Th,mf);
    M = (sqrt(h*Pf*knim*Af))*((mean([Ti,To])- Tb));
    m = sqrt((h*Pf)/(knim*Af));
    Qf = M*tanh(m*hxh/2)*(N-1)
    [QTE P] = Unicouple(Tb,TempC);
    QTET = F*QTE
     if Qf<QTE
            fprintf('Tbg too small\n');
     else
        end
    i=i+1;

end

To = Ti - (Q./(mf*cp));
Th-Tb
    
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
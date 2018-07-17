% Cold Side Heat Exchanger
% Composed of Micro-Plate Fins
% Overall Convection Coefficient obtained from Kim, Park & Lee's paper

clc
clear
tic

hxwC = 0.04;
hxlC = 0.04;
hxhC = 200*(10^-6);
tf_C = 20*(10^-6);

Tb = 30;
Tamb =20;

Rth_TE = 3.901;%2.52; %3.901

NCv = 300:100:2000;
tf_Cv = (0.1*10^-6):(0.1*10^-6):(1*10^-6);
%%
for i = 1:length(NCv)
    for j = 1:length(tf_Cv)
        
    NC = NCv(i);
    tf_C = tf_Cv(j);
    sC = (hxwC - (NC.*tf_C))./(NC-1);  % Fin Spacing
    scm(i,j) = sC;
    
    if sC <= (100*10^-9)
        fprintf('Error!')
        return
    else
    end

    Ar = (hxwC*hxhC*2*NC) + (tf_C*hxwC*NC) + (hxhC*tf_C*NC*2) + sC*hxwC*(NC-1);
    Arm(i,j) =Ar;

    h = 3.99250095E+14*(sC)^4 - 1.92719077E+11*(sC)^3 - 3.37531698E+07*(sC)^2 + 3.67168889E+04*sC + 3.42357891E+00;
    hm(i,j) = h;
    
    Q = h*Ar*(Tb-Tamb);
    Qm(i,j) = Q;
    
    Rth = (Tb-Tamb)/Q;
    Rthm(i,j) =Rth;
    
    Therm_Ratio = Rth/(Rth_TE+Rth);
    Therm_Ratiom(i,j) = Therm_Ratio;
    
    end
end

toc
Min_ThermRes = min(min(Rthm))
Max_h = max(max(hm))
Max_Ar = max(max(Arm))
Min_ThermResRatio = min(min(Therm_Ratiom))
Min_Spacing = min(min(scm))
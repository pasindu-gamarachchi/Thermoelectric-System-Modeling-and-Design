% Overall HT Coeff obtained from [1] Micheli, L., Reddy, K., and Mallick, T., 2015, 
..."General correlations among geometry, orientation and thermal performance of natural
...convective micro-finned heat sinks", International Journal of Heat and Mass Transfer, 91, pp. 711-724.
% Dimensions from paper used
% Date Developed: 06/07


clc
clear

tf_C = 200*(10^-6);
fh = 600*(10^-6);
NCx = 124;
NCy = 124;
hxlC = 0.05;
hxwC = 0.05;

sC = (hxwC - (NCx.*tf_C))./(NCx+1);  % Theoretical Fin Spacing
asC = (hxwC - (NCx.*tf_C))./(NCx-1);  % Actual Fin Spacing


Tb = 30;
Tamb = 20;
x = Tb -Tamb; % Delta T
 

h = 1.84741296E-06*x^3 - 6.18101214E-04*x^2 + 8.54527820E-02*x + 2.66036240E+00;

Ar = (tf_C*fh*4*NCx*NCy) + (tf_C*tf_C*NCx*NCy) + (asC*hxlC*(NCx-1)) +(asC*hxwC*(NCy-1));

Q = Ar*h*(Tb-Tamb);

Rth = x/Q
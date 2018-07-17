% Finite Element for Fin

clc
clear all
close all

% Input
Tb = 30; % Base Temperature
Tamb = 22; % Amb Temp

FinType = 'Pin'
fd = (100*10^-6); % Fin Diameter
fh = (10000*10^-6); % Fin Height
fw = 0.002; % Fin Width
ft = 0.003; % Fin Thickness

TempData = xlsread('Tempr_Midpoint_1213.xlsx','sheet1');
y = TempData(:,1)./1000;
AmbTemp = TempData(:,2);


elemsize = y(2)-y(1);
% elemsize = 10.*(y(2)-y(1));

elems = ceil(fh/elemsize) +1; % Number of Elements
numn =elems+1;


h = 167.65; % Conv Coefficient
k = 396.78; % Thermal Conductivity

l = fh/elems; % Element length
Tambv = AmbTemp(1:numn);
Tamb = mean(Tambv)

% Element K Matrix

if FinType == 'Pin'
    A = pi*(fd^2)*0.25;
    P = pi*fd;
elseif FinType == 'Plate'
    A = fw*ft;
    P = ft*2 + fw*2;
end
% P = 2.8;
C1 = k*A/l;
C2 = h*P*l/6;
C3 = h*P*l*Tamb/2;

ke = C1*[1, -1; -1,1] + C2*[2,1; 1,2];
fe = [C3,C3]';

K = zeros(numn, numn);
F = zeros(numn,1);

for i = 1:elems
    
    dof = [ i, i+1];
    K(dof, dof) = ke +  K(dof, dof);
    F(dof) = fe + F(dof);
    
    
end

%  F(end) = 6.4;
F(2) = F(2) - Tb*(ke(2,1));

for i=1:numn
    K(i,1) =0;
    K(1,i) =0;
    F(1) =Tb;
end

K(1,1)= 1;

Tempr = K\F;
Tempr(1) = Tb;


% Analytical Solution

m = sqrt((h*P)/(k*A));
M = sqrt(h*P*k*A)*(Tb-Tamb);

x = 0:l:fh;

Th_Thb =cosh(m*(fh-x))/(cosh(m*fh));
An_Temp = Th_Thb*(Tb-Tamb) +Tamb;

qf = M*tanh(m*fh);

figure (1)
plot(Tempr, 'bo')
hold on
plot(An_Temp, 'rx')

% Heat Flow

for i = 1:elems
    Tavg(i) = (Tempr(i) +Tempr(i+1))*0.5;
    Q_elem(i) = h*P*l*(Tavg(i) - Tamb);
end

Qf_fe = sum(Q_elem)


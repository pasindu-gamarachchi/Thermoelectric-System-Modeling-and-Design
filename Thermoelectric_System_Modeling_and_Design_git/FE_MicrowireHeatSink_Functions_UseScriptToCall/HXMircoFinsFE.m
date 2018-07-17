% Micro-wire Heat Sink to calculate heat transferred from Heat Sink for 
% varied ambient temperature using finite elements for fins
% Data Created : 01/11/2017

clc
clear all
close all
tic

global  cvl Tcg hxwC  Cb NCx NCz fd fh Tb Tamb elemsize

fd = 100*(10^-6);

cvl = 0.04;
Tcg = 0.01;
hxwC = cvl;
Cb = 0.005;
NCx = 10;
NCz = 10;
fh = 0.005;
Tb = 30;
Tamb = 22;


TemprDat = xlsread('FlatPlateTempr.xlsx', 'sheet1');

y = TemprDat(:,1);

x_100_z_100 = TemprDat(:,2);
x_75_z_100 = TemprDat(:,3);
x_50_z_100 = TemprDat(:,4);
x_25_z_100 = TemprDat(:,5);
x_0_z_100 = TemprDat(:,6);

x_100_z_75 = TemprDat(:,7);
x_75_z_75 = TemprDat(:,8);
x_50_z_75 = TemprDat(:,9);
x_25_z_75 = TemprDat(:,10);
x_0_z_75 = TemprDat(:,11);

x_100_z_50 = TemprDat(:,12);
x_75_z_50 = TemprDat(:,13);
x_50_z_50 = TemprDat(:,14);
x_25_z_50 = TemprDat(:,15);
x_0_z_50 = TemprDat(:,16);

x_100_z_25 = TemprDat(:,17);
x_75_z_25 = TemprDat(:,18);
x_50_z_25 = TemprDat(:,19);
x_25_z_25 = TemprDat(:,20);
x_0_z_25 = TemprDat(:,21);

x_100_z_0 = TemprDat(:,22);
x_75_z_0 = TemprDat(:,23);
x_50_z_0 = TemprDat(:,24);
x_25_z_0 = TemprDat(:,25);
x_0_z_0 = TemprDat(:,26);

Nx = NCx/2;
Nz = NCz/2;

elemsize = y(2)-y(1);

TambM = zeros(Nx,Nz,length(y));

% for j = 1:Nz
for i = 1:Nx
    TambM(i,1,:) = TemprDat(:,i+1);
    TambM(i,2,:) = TemprDat(:,i+6);
    TambM(i,3,:) = TemprDat(:,i+11);
    TambM(i,4,:) = TemprDat(:,i+16);
    TambM(i,5,:) = TemprDat(:,i+21);
end
% end

for i = 1:Nx
    for j = 1:Nz
        
        [Qf_fe, TTip] = FEFinV(TambM(i,j,:));
        Qf_feM(i,j) = Qf_fe;
        TtipM(i,j) = TTip;
    end
end


Qf_Tot = sum(sum((Qf_feM)))
%%
px = (hxwC/2).*[1:-0.25:0]
pz = (hxwC/2).*[1:-0.25:0]

figure(1)
surf(px,pz, TtipM)

toc


for i =1:5
    for j =1:5
        figure(2)
        patch(px(i),pz(j),TtipM(i,j))
        hold on
    end
end


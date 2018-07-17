
global hxwC NCx NCz elemsize

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
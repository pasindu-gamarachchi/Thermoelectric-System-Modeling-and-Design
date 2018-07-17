% Micro-wire Heat Sink to calculate heat transferred from Heat Sink for 
% varied ambient temperature using finite elements for fins
% Data Created : 01/11/2017
% Edited and made into a function using ambient temperature from : 02/09/2017

 function QTot = HXMicroFinsFiniteElement(Tb, Tamb)

global  hxlC Tcg hxwC  Cb NCx NCz fd fh 


temprfile = ['Base', num2str(Tb),'Amb', num2str(Tamb), 'FluidTempDat.mat'];

load(temprfile);


Nx = floor(NCx/2); % Number of fins for a quarter model
Nz = floor(NCz/2); % Number of fins for a quarter model

elemsize = y(2)-y(1); %  added yesterday
elems = ceil(fh/elemsize) -1; % Number of Elements   added yesterday
numn =elems+1; % added yesterday
len = size(Tempr);


Nx_norm = round((1:1:Nx).*(len(1))./(Nx));
Nz_norm = round((1:1:Nz).*(len(2))./(Nz));


%%
TambM = zeros(Nx,Nz,length(y));

% for j = 1:Nz
for i = 1:Nx
    for j = 1:Nz
    
    p = Nx_norm(i);
    q = Nz_norm(j);
    TambM(i,j,:) = Tempr(p,q,:);
   
    end
end
% end
% TemprM = zeros(length(Nx) , length(Nz),  length(numn));
% Tstore  = zeros(length(numn));

for i = 1:Nx
    for j = 1:Nz
        
        [Qf_fe, TTip, Tempr] = FEFinV(TambM(i,j,:),Tb, Tamb, elemsize );
%         for k = 1:length(Tempr)
%             
%             Tstore(k) = Tempr(k);
%         end
%        
        %Tstore  = Tempr
%          TemprM(i,j, :) = Tstore;
        Qf_feM(i,j) = Qf_fe;
        TtipM(i,j) = TTip;
    end
end

% Heat Transfer from Base


h_b = convec_hotplate(Tb, Tamb);
Fin_Ar = 0.25*pi*(fd^2)*(NCx*NCz);
Un_Ar =  (hxwC*hxlC) - Fin_Ar;


Q_b = Un_Ar*(h_b)*(Tb-Tamb); % Heat transfer from Base Area

Qf_Tot = sum(sum((Qf_feM)));
QTot = Qf_Tot*4 + Q_b;


 end



%%
% px = (hxwC/2).*[1:-0.25:0]
% pz = (hxwC/2).*[1:-0.25:0]

% figure(1)
% surf(x,z, TtipM)

% load('Base29.5Amb22FluidTempDat.mat')
% for i =1:5
%     for j =1:5
%         figure(2)
%         patch(px(i),pz(j),TtipM(i,j))
%         hold on
%     end
% end


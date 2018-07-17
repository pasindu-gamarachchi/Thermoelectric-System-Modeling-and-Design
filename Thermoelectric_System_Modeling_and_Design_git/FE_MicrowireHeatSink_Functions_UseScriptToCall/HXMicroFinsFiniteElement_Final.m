% Micro-wire Heat Sink to calculate heat transferred from Heat Sink for 
% varied ambient temperature using finite elements for fins
% Pasindu Gamarachchi - Email : pgamarachchi@gmail.com


 function QTot = HXMicroFinsFiniteElement(Tb, Tamb)

global  hxlC  hxwC   NCx NCz fd fh 


temprfile = ['Base', num2str(Tb),'Amb', num2str(Tamb), 'FluidTempDat.mat']; % Load Temperature File based on base and ambient temperature

load(temprfile);

Nx = floor(NCx/2); 
Nz = floor(NCz/2); 

elemsize = y(2)-y(1);
elems = ceil(fh/elemsize) -1;
numn =elems+1; 
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


for i = 1:Nx
    for j = 1:Nz
        
        [Qf_fe, TTip, Tempr] = FEFinV(TambM(i,j,:),Tb, Tamb, elemsize ); % Function for heat transfer from fins
        Qf_feM(i,j) = Qf_fe;
        TtipM(i,j) = TTip;
    end
end

% Heat Transfer from Base

h_b = convec_hotplate(Tb, Tamb); % Function to obtain convection coeff for Heat Transfer from Base
Fin_Ar = 0.25*pi*(fd^2)*(NCx*NCz);
Un_Ar =  (hxwC*hxlC) - Fin_Ar;


Q_b = Un_Ar*(h_b)*(Tb-Tamb); % Heat transfer from Base Area

Qf_Tot = sum(sum((Qf_feM)));
QTot = Qf_Tot*4 + Q_b;


 end



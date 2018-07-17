% Finite Element for Fin function to be used with heat sink


 function [Qf_fe, Ttip, Tempr] =  FEFinVFTemp(AmbTemp, Tb, Tamb, elemsize)

global fd fh 



elems = ceil(fh/elemsize) -1; % Number of Elements

numn =elems+1;


h = convectionmicro(Tb,Tamb);
k = 396.78; % Thermal Conductivity

l = fh/elems; % Element length
Tambv = AmbTemp(1:numn);

% Element K Matrix

A = pi*(fd^2)*0.25;
P = pi*fd;

C1 = k*A/l;
C2 = h*P*l/6;


ke = C1*[1, -1; -1,1] + C2*[2,1; 1,2];


K = sparse(numn, numn);
F = sparse(numn,1);

for i = 1:elems
    
    dof = [ i, i+1];
    K(dof, dof) = ke +  K(dof, dof);
 
   
end

% Calculates Thermal Load for each node
for i = 1:numn
    fe(i) = h*P*l*Tambv(i)/2;
    
end


% Combines the Thermal Loads from previous element
for i = 2:numn-1
    F(i) = fe(i)+ fe(i-1);
    
end

F(end) = fe(end); % Last node only has thermal load from itself
F(2) = F(2) - Tb*(ke(2,1));
 
for i=1:numn
    K(i,1) =0;
    K(1,i) =0;
    F(1) = Tb;
end

K(1,1) =1;

Tempr = K\F;


for i = 1:elems
    Tavg(i) = (Tempr(i) +Tempr(i+1))*0.5;
    Q_elem(i) = h*P*l*(Tavg(i) - Tambv(i));
end


Qf_fe = sum(Q_elem);
Ttip = Tempr(end);


 end

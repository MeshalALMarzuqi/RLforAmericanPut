%load('C:\Analysis\Data\Results_0');

%Obj             = Real_Payoffs;
Obj             = Sim_Payoffs;
Size            = size(Obj);
I               = Size(4);
Eff_Sample_Size = zeros(Size(1),Size(2),Size(3));

for k1 = 1:Size(1)
    for k2 = 1:Size(2)
        for k3 = 1:Size(3)
            Eff_Sample_Size(k1,k2,k3)   = sum(~isnan(Obj(k1,k2,k3,:)));   
        end 
    end
end
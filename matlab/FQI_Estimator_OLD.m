function [Weights,Converged,Iter] = FQI_Estimator_OLD(Ini_Weights,Data,Vol_Data,r,Settings)
Use_Vol     = (size(Vol_Data,1) ~= 0);
W           = 7 + 2 * Use_Vol;              % The number of weights depends on whether we use volatility or not.
if size(Ini_Weights,1) > 0
    Weights     = Ini_Weights;    
else
    Weights     = zeros(W,1);
end
gamma       = exp(-r * Settings.dt);
Tenor       = size(Data,2) - 1;
N           = size(Data,1);
Basis       = zeros(W,N,1+Tenor);
A           = zeros(W,W);

% Calculate A offline, since it does not change with Weights.
for j = 1 : N
    for t = 0 : (Tenor-1)
        S_t         = Data(j,1+t);
        Phi_t       = [1; exp(-S_t/2); exp(-S_t/2)*(1-S_t); exp(-S_t/2)*(1-2*S_t+S_t^2/2); sin(-t/Tenor * pi/2 + pi/2); log(Tenor - t); (t/Tenor)^2];            
        Basis(:,j,1+t)  = Phi_t;
        A               = A + Phi_t * Phi_t';  % W-by-W matrix.
    end
end
Inv_A       = inv(A);

% Calculate b online, since it is affected by Weights.
Converged   = 0;
Iter        = 1;
while (Iter < max(Settings.Min_Iter,10) || ~Converged) && Iter <= Settings.Max_Iter
    Old_Weights = Weights;
    b           = zeros(W,1);        
    for j = 1 : N_Mod
        for t = 0 : (Tenor-1)
            S_tp1       = Data(j,1+t+1);
            Phi_t       = Basis(:,j,1+t);
            Phi_tp1     = Basis(:,j,1+t+1);            
            g_tp1       = max(1 - S_tp1, 0); % The exercise value of a put option with strike 1. 
            % We use strike 1 because the initial price is 1, so 1 is the ATM strike.
            b           = b + Phi_t * max(g_tp1, Phi_tp1' * Weights); % W-by-1 vector.
        end
    end
    if abs(det(A)) > Settings.Det_Tol^W
        Weights     = Inv_A * b;
    end 
    Converged   = max( 2 * abs(Weights - Old_Weights) ./ (Weights + Old_Weights) ) <= Settings.Weight_Tol;
    Iter        = Iter + 1;  
end
Iter        = Iter - 1;
if ~Converged && strcmpi(Settings.Warnings,'on')
    fprintf( '\n                  FQI algorithm has not converged.' );
end
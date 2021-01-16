%% Longstaff-Schwartz Method
function [Weights,price] = LSM_Estimator(Data,r,Settings)
price       = 0;
W           = 4;              % The number of weights depends on whether we use volatility or not.
Tenor       = size(Data,2) - 1;
N           = size(Data,1);
Weights     = zeros(W,1+Tenor);
gamma       = exp(-r * Settings.dt);
Converged   = 0;
Iter        = 1;
while (Iter < Settings.Min_Iter || ~Converged) && Iter <= Settings.Max_Iter
    Old_Weights = Weights;
    A           = zeros(W,W);
    b           = zeros(W,1);    
    for t = (Tenor-1) : -1 : 0
        for j = 1 : N     
            S_t         = Data(j,1+t);
            S_tp1       = Data(j,1+t+1);
            Phi_t       = [1; exp(-S_t/2); (exp(-S_t/2))*(1-S_t); (exp(-S_t/2))*(1 - 2*S_t + S_t^2/2)];   
            if t == Tenor - 1
                Phi_tp1     = zeros(W,1);
            else
                Phi_tp1     = [1; exp(-S_tp1/2); exp(-S_tp1/2)*(1-S_tp1); exp(-S_tp1/2)*(1-2*S_tp1+S_tp1^2/2)];               
            end            
            g_tp1 = max(1 - S_tp1, 0); % The exercise value of a put option with strike 1. 
            % We use strike 1 because the initial price is 1, so 1 is the ATM strike.
            A     = A + Phi_t * Phi_t';                                          % W-by-W matrix.
            b     = b + gamma * Phi_t * max(g_tp1, Phi_tp1' * Weights(:,1+t+1)); % W-by-1 vector.
        end       
        if abs(det(A)) > Settings.Det_Tol^W
            Weights(:,1+t)  = inv(A) * b;
        end
    end   
    Converged   = max(max( 2 * abs(Weights - Old_Weights) ./ (Weights + Old_Weights) )) <= Settings.Weight_Tol;
    Iter        = Iter + 1;
end
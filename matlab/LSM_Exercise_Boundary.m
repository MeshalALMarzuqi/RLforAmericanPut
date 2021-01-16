function [Boundary] = LSM_Exercise_Boundary(Weights,Data,sigma,r,Settings)            

% Compiling the exercise boundary for the first path. The boundary B_t is solved from the following inequality: 
%
%   max(1 - S_t, 0) = Intrinsic Value > Intrinsic Value Boundary = max(1 - B_t, 0) = Continuation Value(B_t) = Phi_t(B_t)' * Weights
%

Use_Vol     = size(sigma,1) == 1 && size(sigma,2) == 1 && ~isnan(sigma);
Tenor       = size(Data,2) - 1;
Boundary    = zeros(1,Tenor) * NaN;
S_Grid      = 0.1:0.001:3;
I           = length(S_Grid);
ind         = 1:I;

for t = 0:(Tenor-1)
    I_Vals          = zeros(I,1) * NaN;
    C_Vals          = zeros(I,1) * NaN;
    
    for i = 1:I
        S_t         = S_Grid(i);
        Intr_Value  = max(1 - S_t, 0);                  % The intrinsic (exercise) value of a put option with strike 1. 
        
        Phi_t       = [1; exp(-S_t/2); exp(-S_t/2)*(1-S_t); exp(-S_t/2)*(1-2*S_t+S_t^2/2)];
        if Use_Vol && Settings.Vol_Basis == 1           % Use simple basis functions based on volatility.
            Phi_t       = [Phi_t; sigma * sqrt(Tenor - t) * S_t; sigma^2 * (Tenor - t) * S_t];
        elseif Use_Vol                                  % Use basis functions which are based on volatility and inspired by linearization of the Black-Scholes formula.
            sigma_sq    = sigma^2;
            tau         = Tenor - t;
            tau_sqrt    = sqrt(tau);
            d1          = (log(S_t) + (r + sigma_sq / 2) * tau) / (sigma * tau_sqrt);
            NPrime      = exp(- d1^2 / 2);
            Phi_t       = [Phi_t; sigma * tau_sqrt * S_t * NPrime; sigma_sq * tau * S_t * NPrime];
        end         
        
        Cont_Value  = Phi_t' * Weights(:,1+t);
        I_Vals(i)   = Intr_Value;
        C_Vals(i)   = Intr_Value - Cont_Value; % our f(S) equation
    end
    
    flags           = I_Vals > 0 & C_Vals > 0;
    if sum(flags) > 0
        Boundary(t+1)   = S_Grid(max(ind(flags)));
    end
end
% figure, plot(Boundary); xlabel('Time[years]'); ylabel('Stock Price [a.u.]'); title('LSM, Exercise Boundary, Exhaustive')
end
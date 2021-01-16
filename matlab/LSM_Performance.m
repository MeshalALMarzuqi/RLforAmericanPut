function [Payoff_Mean,Payoff,Cont_Boundary] = LSM_Performance(Weights,Data,Vol_Data,r,Settings, Ticker, MARKER_plot, method, K)
tol           = 0.001;
Use_Vol     = (size(Vol_Data,1) ~= 0);
Tenor       = size(Data,2) - 1;
N           = size(Data,1);
Payoff      = zeros(N,1) * NaN;
Cont_Boundary = zeros(Tenor,1) * NaN;

% Trading begins.
for j = 1:N
    Exercised   = 0;
    t           = 0;
    
    while ~Exercised && t <= Tenor - 1
        S_t         = Data(j,1+t);
        Intr_Value  = max(1 - S_t, 0);                  % The intrinsic (exercise) value of a put option with strike 1.
        
        Phi_t       = [1; exp(-S_t/2); exp(-S_t/2)*(1-S_t); exp(-S_t/2)*(1-2*S_t+S_t^2/2)];
        if Use_Vol && Settings.Vol_Basis == 1           % Use simple basis functions based on volatility.
            Phi_t       = [Phi_t; Vol_Data(j,1+t) * sqrt(Tenor - t) * S_t; Vol_Data(j,1+t)^2 * (Tenor - t) * S_t];
        elseif Use_Vol                                  % Use basis functions which are based on volatility and inspired by linearization of the Black-Scholes formula.
            sigma       = Vol_Data(j,1+t);
            sigma_sq    = sigma^2;
            tau         = Tenor - t;
            tau_sqrt    = sqrt(tau);
            d1          = (log(S_t) + (r + sigma_sq / 2) * tau) / (sigma * tau_sqrt);
            NPrime      = exp(- d1^2 / 2);
            Phi_t       = [Phi_t; sigma * tau_sqrt * S_t * NPrime; sigma_sq * tau * S_t * NPrime];
        end
        
        Cont_Value  = Phi_t' * Weights(:,1+t);
        Exercised   = (Intr_Value > Cont_Value);
        
        if Exercised
            Payoff(j)   = exp(-r * t * Settings.dt) * Intr_Value;                   % Discounted payoff.
        else
            t           = t + 1;
        end
    end
    
    if ~Exercised
        Payoff(j)   = exp(-r * Tenor * Settings.dt) * max(1 - Data(j,1+Tenor), 0);  % Discounted payoff.
    end
end
Payoff_Mean     = mean(Payoff);

%% CALCULATING THE BOUNDARY
% cont. boundary
for t = 1 : Tenor
    midS = 0;
    midf = 0;
    for midS = 1 : -0.0001 : 0 % have to go from 1 down, since need to find the 1st stock price for which the midf=0. Since all stocks below it will give f(S) = 0 also!
        %------------ calc Q(S)
        Phi_t       = [1; exp(-midS/2); exp(-midS/2)*(1-midS); exp(-midS/2)*(1-2*midS+midS^2/2)];
        Cont_Value  = Phi_t' * Weights(:,t); % Q prediction
        %------------ end calc Q(S)
        midf = fBoundary(midS,Cont_Value,1); % midf = S - K + Q; NOTE: Cont_Value is a function of midS!
        if abs(midf) < tol
            Cont_Boundary(t,1) = midS*K;
            break
        end
    end % end while
end   
Cont_Boundary(end,1) = K;
% [Cont_Boundary_interp] = Data_Interpolation(Cont_Boundary);

if MARKER_plot
    figure, plot(Cont_Boundary, 'x-', 'MarkerSize', 18, 'LineWidth', 4); xlabel('Time [days]'); ylabel('Stock Price [a.u.]'); 
    title(['LSM. Company: ' Ticker, ' Tenor ,' num2str(Tenor), ' Method: ' method])
%     OutputName = Ticker; % Tikcer is already a string, hence don't need ''
%     saveas(h, ['figures/' OutputName '.fig']); 
%     OutputName = Ticker;
%     save(['synth1/memo/' OutputName '.mat'], 'Var1', 'Var2');
%     saveas(FigureNumber, ['synth1/memo/' OutputName '.jpg']);
end    
end

function [sol] = fBoundary(S,Q,K)
sol = S-K+Q;
end
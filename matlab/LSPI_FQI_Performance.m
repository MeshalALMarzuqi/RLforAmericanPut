function [Payoff_Mean,Payoff,Cont_Boundary] = LSPI_FQI_Performance(Weights,Data,Vol_Data,r,Settings, Ticker, MARKER_plot, method, K)

% Cont_Value is Q(t,S,sigma) is the boundary we believe.
tol = 0.001;
Tenor         = size(Data,2);
N             = size(Data,1);
Payoff        = zeros(N,1) * NaN;
Cont_Boundary = zeros(Tenor,1) * NaN;

for j = 1 : N % number of trajs of stock rollouts
    Exercised = 0;
    t         = 2;   
    while ~Exercised && t <= Tenor
        S_t           = Data(j,t);
        Intr_Value    = 1 - S_t;                               % The intrinsic value of a put option with strike 1.
        if Intr_Value > 0
            Phi_t         = [1; exp(-S_t/2); exp(-S_t/2)*(1-S_t); exp(-S_t/2)*(1-2*S_t+S_t^2/2); sin(-t/Tenor * pi/2 + pi/2); log(Tenor - t); (t/Tenor)^2];
            Exercised     = Intr_Value > (Phi_t' * Weights);   % K - S(t) > Q(S,t)     
        end
        if Exercised
            Payoff(j)     = exp(-r * t * Settings.dt) * Intr_Value;                % Payoff discounted to present time
        else
            t             = t + 1;
        end
    end   
    if ~Exercised
        Payoff(j)         = exp(-r*Tenor*Settings.dt) * max(1 - Data(j,Tenor), 0); % Payoff discounted to present time
    end
end

%% FOR CALCULATING THE BOUNDARY
% --- EXHAUSTIVE root searching algorithm!!!
% search all stock prices from 1 down to 0 and find 1st stock price for which f(S) = S - K + Q = 0
for t = 1 : Tenor
    midS = 0;
    midf = 0;
    for midS = 1 : -0.0001 : 0 % have to go from 1 down, since need to find the 1st stock price for which the midf=0. Since all stocks below it will give f(S) = 0 also!
        %------------ calc Q(S)
        Phi_t       = [1; exp(-midS/2); exp(-midS/2)*(1-midS); exp(-midS/2)*(1-2*midS+midS^2/2); sin(-t/Tenor * pi/2 + pi/2); log(Tenor - t); (t/Tenor)^2];
        Cont_Value  = Phi_t' * Weights; % Q prediction
        %------------ end calc Q(S)
        midf = fBoundary(midS,Cont_Value,1); % midf = S - K + Q; NOTE: Cont_Value is a function of midS!
        if abs(midf) < tol
            Cont_Boundary(t,1) = midS*K;
            break
        end
    end % end while
    
% --- BISECTION root searching algorithm!!!
% %     lS = 0.60;
% %     uS = 0.95;
% %     tol = 10^(-12);
% %     %lf = fBoundary(lS,Cont_Value,1); % expect to be < 0
% %     %uf = fBoundary(uS,Cont_Value,1); % expect to be > 0
% %     midS = 0;
% %     midf = 0;
% %     while uS - lS > tol
% %         midS = (uS + lS) * 0.5;
% %         %------------ calc Q(S)
% %         Phi_t       = [1; exp(-midS/2); exp(-midS/2)*(1-midS); exp(-midS/2)*(1-2*midS+midS^2/2); sin(-t/Tenor * pi/2 + pi/2); log(Tenor - t); (t/Tenor)^2];
% %         Cont_Value  = Phi_t' * Weights; % Q prediction 
% %         %------------ end calc Q(S)
% %         midf = fBoundary(midS,Cont_Value,1); % Cont_Value is a function of midS
% %         if midf > 0
% %             uS = midS;
% %             %uf = midf;
% %         else
% %             lS = midS;
% %             %lf = midf;
% %         end
% %     end % end while
% %     if midS > 0.61
% %         Cont_Boundary(t,1) = midS;
% %     end
end
Cont_Boundary(end,1) = K;
if MARKER_plot
    figure, plot(Cont_Boundary, 'x-', 'MarkerSize', 18, 'LineWidth', 4); xlabel('Time [days]'); ylabel('Stock Price [a.u.]'); axis tight;
    title(['LSPI. Company: ' Ticker, ' Tenor ,' num2str(Tenor), ' Method: ' method])
%     OutputName = Ticker; % Tikcer is already a string, hence don't need ''
%     saveas(h, ['figures/' OutputName '.fig']); 
%     OutputName = Ticker;
%     save(['synth1/memo/' OutputName '.mat'], 'Var1', 'Var2');
%     saveas(FigureNumber, ['synth1/memo/' OutputName '.jpg']);
end
Payoff_Mean    = mean(Payoff); % average payoff of all trajs
end

function [sol] = fBoundary(S,Q,K)
sol = S - K + Q; % i.e. K - S - Q > 0 is exercise. Thus continuation is when we keep i.e. - K + S + Q < 0.
end
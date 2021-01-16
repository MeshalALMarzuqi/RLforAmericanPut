function [Payoff_Mean,Payoff,Cont_Boundary_interp] = LSM_Performance_OriginalPaper(Weights,Data,Vol_Data,r,Settings, Ticker, MARKER_plot, method, K)
numBasis_S    = 4;
[N,Tenor]     = size(Data); 
tol           = 0.001;
Payoff        = zeros(N,1);
Cont_Boundary = zeros(Tenor,1) * NaN;

%% LSM TEST
for ii = 1 : N % num trajs
    Phi_lsm = zeros(numBasis_S, Tenor); % ### THEY DONT USE TIME BASIS FOR LSM FOR SOME REASON
    trajectory_test = Data(ii,1:Tenor-1);    
    Phi_lsm(:,1:Tenor-1) = [repmat(1,Tenor-1,1)'; exp(-trajectory_test/2); exp(-trajectory_test/2).*(1-trajectory_test); exp(-trajectory_test/2).*(1-2*trajectory_test+trajectory_test.^2/2)];
%     for k = 1:numBasis_S
%         Phi_lsm(k,1:Tenor-1) = feval(bases{k}, trajectory_test(1:Tenor-1));
%     end
    for t = 2 : Tenor
        intrinsicValue = 1 - Data(ii,t);
        if intrinsicValue > 0
            continuationValue = Phi_lsm(:,t)' * Weights(:,t);
            if intrinsicValue > continuationValue
                Payoff(ii,1) = intrinsicValue * exp(-r*Settings.dt*t);
                break
            end
        end % end intrinsicValue
    end % end Tenor
%     if isnan(Payoff(ii,1)) % if not exercised prior to expiry, value is like european.
%         Payoff(ii,1) = exp(-r * Tenor * Settings.dt) * max(1 - Data(ii,Tenor), 0);  % Discounted payoff.
%     end
end
Payoff_Mean = mean(Payoff);
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
Cont_Boundary(53,1) = 22;
Cont_Boundary(60,1) = 20;
Cont_Boundary(70,1) = 18;
Cont_Boundary(end,1) = K;
[Cont_Boundary_interp] = Data_Interpolation(Cont_Boundary);

if MARKER_plot
    figure, plot(Cont_Boundary_interp, 'x-', 'MarkerSize', 18, 'LineWidth', 4); xlabel('Time [days]'); ylabel('Stock Price [a.u.]'); 
    title(['LSM. Company: ' Ticker, ' Tenor ,' num2str(Tenor), ' Method: ' method])
%     OutputName = Ticker; % Tikcer is already a string, hence don't need ''
%     saveas(h, ['figures/' OutputName '.fig']); 
%     OutputName = Ticker;
%     save(['synth1/memo/' OutputName '.mat'], 'Var1', 'Var2');
%     saveas(FigureNumber, ['synth1/memo/' OutputName '.jpg']);
end     
end

function [sol] = fBoundary(S,Q,K)
sol = S - K + Q; % i.e. K - S - Q > 0 is exercise. Thus continuation is when we keep i.e. - K + S + Q < 0.
end
function [Weights,price] = FQI_Estimator(Ini_Weights,Data,Vol_Data,r,Settings)
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
last_price = -99; % initialise var used to determine option price convergence
for iter = 1 : Settings.Max_Iter
    b = zeros(W,1);
    for j = 1 : N
        for t = 0 : (Tenor-1)
            S_tp1   = Data(j,1+t+1);
            Phi_t   = Basis(:,j,1+t);
            Phi_tp1 = Basis(:,j,1+t+1);
            g_tp1   = max(1 - S_tp1, 0); % The exercise value of a put option with strike 1.
            % We use strike 1 because the initial price is 1, so 1 is the ATM strike.
            b       = b + Phi_t * max(g_tp1, Phi_tp1' * Weights); % W-by-1 vector.
        end
    end
    if(rank(A) == W)
        Weights = A\b;
    else
        Weights = pinv(A)*b;
    end
    % ----- Calculate option price after 1 iteration of FQI algorithm
    payoff = zeros(N,1);
    for ii = 1 : N;
        for t = 2 : Tenor
            intrVal = 1 - Data(ii,t); % standardised strike - S_t
            if intrVal > 0
                Q = Basis(:,ii,t)'*Weights;
                if intrVal > Q
                    payoff(ii) = exp(-r*t*Settings.dt) * intrVal; % discount intrVal to present time
                    break % if intrVal > Q execute. Terminate current for loop (current traj loop) --> go onto evaluating next traj
                end
            end
        end
    end
    price = mean(payoff);
    price_Diff = abs(price - last_price); % MEAN PAYOFF + 99!!!!
    if (iter > 3)  &&  (price_Diff < Settings.price_Tol) % perform at least 3 iterations
        fprintf(['Converged at Iteration: ' num2str(iter) '\n'])
        break
    end
    last_price = price;
    if iter == Settings.Max_Iter
        fprintf(['LSPI ran for MAX # Iterations: ' num2str(iter) '\n'])
    end
end
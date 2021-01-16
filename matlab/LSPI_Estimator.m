%% Least Squares Policy Iteration (LSPI) 
% Basis functions f(S,t) depend on S and t, but weights w(t) = w are constant

function [Weights,price] = LSPI_Estimator(Ini_Weights,Data,Vol_Data,r,Settings)
% Vol_Data is the historical (realised) volatility rescaled by dt
Use_Vol = (size(Vol_Data,1) ~= 0); % if Vol_Data is empty (controlled by Vol_Window), then we are not using volatility. True if size(Vol_Data,1) does not equal 0 --> give Use_Vol = 1.
W = 7 + 2 * Use_Vol; % The number of weights depends on whether we use volatility or not: W=9 if use vol; W=7 if not use vol.
if size(Ini_Weights,1) > 0
    Weights = Ini_Weights;    
else
    Weights = zeros(W,1);
end
gamma = exp(-r * Settings.dt); % continuous discount rate
Tenor = size(Data,2); % # time steps in 1 traj
N = size(Data,1); % number of windowed trajs obtaiend from 1 real path
Basis = zeros(W,N,Tenor);
for ii = 1 : N % looping over each trajectory (total N trajs)
    for t = 1 : Tenor-1 % looping over individual trajectories
        S_t = Data(ii,t); % Data is rows (individual traj), columns is time steps
        % feeding states through 7 no vol basis (action not fed through here)
        Basis(:,ii,t) = [1; exp(-S_t/2); (exp(-S_t/2))*(1-S_t); (exp(-S_t/2))*(1-2*S_t+S_t^2/2);... % underlying basis
                        sin(-(t*pi)/(2*Tenor) + pi/2); log(Tenor - t); (t/Tenor)^2];                % time basis
    end
end
% OutputName = 'OurCodeBasis';
% save(['results/' OutputName '.mat'], 'Basis');

last_price = -99; % initialise var used to determine option price convergence
for iter = 1 : Settings.Max_Iter
    A = zeros(W,W);
    b = zeros(W,1);
    for i = 1 : N
        S = Data(i,:)';
        Phi_curr(:, :)         = Basis(:,i,:);
        currQ                  = Phi_curr'*Weights;      
        Phi_next               = Phi_curr;
        Phi_next(:, 1:Tenor-1) = Phi_curr(:, 2:Tenor);
        Phi_next(:, Tenor)     = 0;      
        reward                 = zeros(Tenor, 1);
        % stopIndex              = find( 1-S > currQ ); % stopIndex acts as an actions, since it is later used to make Phi_next = 0 if immidiate payoff (K-S) is > than continuation (currQ) then exercise and need to use 0!
        stopIndex              = find( max(1-S,0) > currQ ); % only exercise if in-the-money
        reward(stopIndex, 1)   = max( 0,1-S(stopIndex) );
        Phi_next(:, stopIndex) = 0; %if K-S > currQ, MEANS EXERCISE. Whether in or out of the money Phi_next doesn't exist. Note the reward at exercise is the intrinsic value. So it could be 0 if out of money or +ve if in the money
        %
        A = A + Phi_curr * (Phi_curr- gamma * Phi_next)';
        %
        b = b + Phi_curr* reward;
    end    
% %     fprintf(['Iteration # ' num2str(iter) '\n'])
% %     A = zeros(W,W);
% %     b = zeros(W,1);
% %     for i = 1 : N % going over all trajectories
% %         % Incremental update rule of LSTDQ to update A,b at each t
% %         for t = 1 : Tenor % going over each traj 1 time step at a time
% %             Phi_t = Basis(:,i,t); % Basis at t
% %             if t == Tenor % DO NOT update b at expiry (t==Tenor) (AS PER PAPER'S CODE)
% %                 Phi_tp1 = zeros(W,1); % day after expiry no Phi (Phi_tp1 IS 0 AT t=Tenor-1 and at t=Tenor, AS PER PAPER'S CODE)
% %             else
% %                 Phi_tp1 = Basis(:,i,t+1); % Basis at t+1 
% %             end
% %             S_t = Data(i,t); 
% %             % ----- The exercise value of a put option with strike 1. We use strike 1 because the initial price is 1, so 1 is the ATM strike.            
% %             % g_tp1       = max(1 - S_tp1, 0); % immidiate reward (not S_tp because it is the reward obtained for going from s to s')                                 
% %             g_t = max(1 - S_t, 0);
% %             % ----- (Policy Pi) Deterministic exercise strategy (used as the action instead of feeding action as Basis(s,a))
% %             % policy improved via Weights upon each itertation of LSPI; Pi True =1; False =0; % currQ = (Phi_t' * Weights);         
% %             Pi_t = 1 - S_t > (Phi_t' * Weights); % Note will exercise if IntrVal > Q (but can have situations where 0 > -ve Q; will still exercise but reward g_t=0)
% %             A = A + Phi_t * (Phi_t - gamma * (1 - Pi_t) * Phi_tp1)'; % W-by-W matrix
% %             % b is 0 during continuation and payoff at exercise. On expiry B should not be updated according to paper's code.
% %             b = b + Phi_t * Pi_t * g_t; % W-by-1 vector. Basis is timed by indicator function
% %         end   
% %     end
    % ----- Update weights after 1 iteration over all trajectories
    if(rank(A) == W)
        Weights = A\b;
    else
        Weights = pinv(A)*b;
    end
    % ----- Calculate option price after 1 iteration of LSPI algorithm
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
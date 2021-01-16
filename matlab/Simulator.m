%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATING SHOCKS FOR THE SAMPLE PATHS

function [Paths] = Simulator(N_Paths,Time_Steps,Model,Model_Param,r,Seed,Settings)
% GBM and GARCH start at S=1 for each path at t_1 (hard coded in Simulator)
% GBM and GARCH are segmented into Tennor lengths upon creation.
if Seed >= 0
    randn('state',Seed);        % Initialize the random number generator.
end
Shocks          = randn(N_Paths,Time_Steps);
Paths           = ones(N_Paths,Time_Steps);

if strcmpi(Model,'GARCH')
    % Initializing the volatility as the long term mean.
    Resc_Vol            = ones(N_Paths,1) * sqrt(Model_Param(1) / (1 - Model_Param(2) - Model_Param(3)));
    for t = 1 : Time_Steps-1
        Paths(:,t+1)    = Paths(:,t) .* exp(r * Settings.dt - Resc_Vol.^2 / 2 + Resc_Vol .* Shocks(:,t));
        Resc_Vol        = sqrt( Model_Param(1) + Model_Param(2) * Resc_Vol.^2 + Model_Param(3) * (Resc_Vol .* Shocks(:,t)).^2 );
    end
else    % Model = 'GBM'
    N = N_Paths/2;
    BS_mu    = ( r - 0.5*Model_Param^2 ) * Settings.dt;
    BS_sigma = Model_Param*sqrt(Settings.dt);
    for t = 1 : Time_Steps-1
        randnum = randn(N, 1);
        % first half of data (usual GBM)
        increment = exp(BS_mu + BS_sigma*randnum);
        Paths(1:N,t+1) = Paths(1:N,t).*increment(:);
        % second half of data has vol taken away! (antithetic variance reduction)
        increment = exp(BS_mu + BS_sigma*(-randnum));
        Paths(N+1:N_Paths,t+1) = Paths(N+1:N_Paths,t).*increment(:);
    end
end
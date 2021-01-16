clearvars
close all
% ---------------------- SELECT ALGORITHM FOR TESTING ---------------------
% m = 1; % LSPI
% m = 2; % FQI
% m = 3; % LSM
% ---------------------- SELECT ALGORITHM FOR TESTING ---------------------

% use either a random or given seed as entry point to random num generator
useSeeding = 1; % will know which seed was used
giveSeed = 1;   % give a particular seed to be used
givenSeed = 112212;
seed = 0;
if useSeeding
    seed = randseed();
    if giveSeed
        seed = givenSeed; % 68903 breaks on 10 subT
    end
    display([ 'seeding experiment with: ' num2str(seed)]);
    s = RandStream('mcg16807','Seed',seed);
    oldStream = RandStream.setGlobalStream(s); % gives out old stream, but replaces it with the specified seed from above
end

r                   = 0.03; % 3% interest rate used for all data
Tenors              = [63, 126, 252]; % maturity of option
N_Training_Paths    = 5000; % MUST BE DIVISIBLE BY 2! Num GBM training paths
N_Test_Paths        = 1000;
Realized_Vol_Window = 1; % 0, 1 or 21 % if <=1 means we are not using volatility. 21 is chosen due to num of business days.
if Realized_Vol_Window == 0 || Realized_Vol_Window == 1
    numBasis = 7; 
elseif Realized_Vol_Window == 21
    numBasis = 9; 
else
    fprintf('Mistake in setting Realized_Vol_Window, has to be 0, 1 or 21')
end


% figure
% for m = 1 : 3
    [Real_Payoffs,Sim_Payoffs,Mu,GBM_Vol,GARCH_Param, results] = American_Option_Pricing(numBasis,Tenors,N_Training_Paths,N_Test_Paths,Realized_Vol_Window,r,m);
%     if m == 1
%         p1 = plot(results.Real_Boundaries{2,13,1},'bx-', 'MarkerSize', 18, 'LineWidth', 4); hold on
%     elseif m == 2
%         p2 = plot(results.Real_Boundaries{2,13,1},'gx-', 'MarkerSize', 18, 'LineWidth', 4); hold on
%     else %lsm
%         p3 = plot(results.Real_Boundaries{2,13,1},'rx-', 'MarkerSize', 18, 'LineWidth', 4); hold on
%     end
% end
% xlabel('Time [days]'); ylabel('Stock Price [a.u.]');
% legend([p1,p2,p3],'LSPI','FQI','LSM')
% axis tight

display([ 'Experement was seeded with: ' num2str(seed)])

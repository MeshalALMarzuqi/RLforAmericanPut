%% Longstaff-Schwartz Method
function [w,price] = LSM_Estimator_Originalpaper(Data,r,Settings)
K              = 1;
numBasis       = 4;             
[N, Tenor]     = size(Data);
w              = zeros(numBasis,Tenor);
discountVector = exp( -r * Settings.dt * (1:Tenor)' );
cashFlows      = max(0, K-Data(:,Tenor)); % 500 x 1
exerciseTime   = Tenor*ones(N, 1);        % start from expiry

for step = (Tenor-1) : -1 : 1
  Q = []; X = []; Y = []; Basis = []; numPts = [];
  inMoney = find( Data(:,step) < K ); % inMoney gives the path numbers which are in the money at time step 'step'
  X       = Data(inMoney,step);       % find all in-the-money paths at time step 'step'
  Y       = cashFlows(inMoney).*discountVector(exerciseTime(inMoney)-step);
  numPts  = length(X);
  Basis   = [ones(numPts,1), exp(-X/2), (exp(-X/2)).*(1-X), (exp(-X/2)).*(1 - 2*X + X.^2/2)]; % Basis [length(X) x numBasis]
  % perform regression 
  if(rank(Basis) == numBasis) 
      w(:,step) = Basis \ Y; 
  else
      w(:,step) = pinv(Basis) * Y; % regression. w_t [numBasis x 1]
  end 
  % update cash flows
  intrinsicValue              = K - X;                    % [length(X) x 1]
  Q                           = Basis * w(:,step);        % [length(X) x 1] i.e. Basis * w; (Continuation value Q); 
  index                       = find(intrinsicValue > Q); % [length(X) x 1]
  exercisePaths               = inMoney(index);           % [length(X) x 1]
  cashFlows(exercisePaths)    = intrinsicValue(index);    % 500 x 1. Overwrite the cashFlows of paths that were exercised. i.e. always keeping the latest info on cashflow
  exerciseTime(exercisePaths) = step;                     % keep track of time step at which exercise happened for discounting later.
end

price = max(0, mean(cashFlows.*discountVector(exerciseTime))); % K-S0 = 0; since pricing at the money option. 
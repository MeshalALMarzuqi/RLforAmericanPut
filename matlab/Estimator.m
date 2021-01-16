function [Weights,price] = Estimator(Ini_Weights,Data,r,Method,Vol_Window,Initial_Vol,Settings)

% Calculating realized volatility rescaled by dt.
Vol_Data = Volatility_Data(Data,Vol_Window,Initial_Vol,Settings);

% Calling a specific method.
if strcmpi(Method,'LSPI')
        [Weights,price] = LSPI_Estimator(Ini_Weights,Data,Vol_Data,r,Settings);
elseif strcmpi(Method,'FQI')
        [Weights,price] = FQI_Estimator(Ini_Weights,Data,Vol_Data,r,Settings);
else    % Method = 'LSM'
        [Weights,price] = LSM_Estimator_Originalpaper(Data,r,Settings);
        % [Weights,price] = LSM_Estimator(Data,r,Settings);
end
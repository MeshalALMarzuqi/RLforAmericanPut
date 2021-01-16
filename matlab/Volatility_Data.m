function [Vol_Data] = Volatility_Data(Data,Vol_Window,Initial_Vol,Settings)

%
% Calculating realized annualized volatility.
%

if Vol_Window <= 1
    Vol_Data            = [];
else
    Tenor               = size(Data,2) - 1;
    Vol_Window_Mod      = min(Vol_Window,Tenor);                % Just in case the tenor of the option is too small.
    Eff_Decay           = log(2) / (Vol_Window_Mod - 1);    
    Lags                = 0:(Vol_Window_Mod-1);
    Var_Weights         = exp(-Eff_Decay * Lags);
    Var_Weights         = Var_Weights / sum(Var_Weights);
    N                   = size(Data,1);
    Vol_Data            = Data * 0;    
    Log_Return          = log(Data(:,2:(1+Tenor))) - log(Data(:,1:Tenor));
    Mean                = mean(mean(Log_Return));
    Sq_Dist             = (Log_Return - Mean).^2;
    
    for t = Vol_Window_Mod:Tenor
        Var_Snapshot        = zeros(N,1);
        
        for s = Lags
            Var_Snapshot    = Var_Snapshot + Var_Weights(1+s) * Sq_Dist(:,t - s);            
        end
        
        Vol_Data(:,1+t) = sqrt(Var_Snapshot) / sqrt(Settings.dt);   % Note that the volatility is annualized.
    end 
    
    for t = 0:(Vol_Window_Mod-1)
        if Initial_Vol < 0
            Vol_Data(:,1+t) = Vol_Data(:,1+Vol_Window_Mod);         % Initially, the volatility is constant for a while.
        else
            Vol_Data(:,1+t) = Initial_Vol;                          % The first several values of volatility are set based on the parameters
                                                                    % obtained during the estimation stage.
        end
    end
end
function [Ave_Payoff,Payoffs,ExBoundary_root,ExBoundary_exh] = Performance(Method,Weights,Data,r,Vol_Window,Initial_Vol,MARKER_Exer_Boundary,Settings, Ticker, MARKER_plot, method, K)
% Calculating realized volatility rescaled by dt.
Vol_Data  = Volatility_Data(Data,Vol_Window,Initial_Vol,Settings);
ExBoundary_exh = [];
if strcmpi(Method,'LSPI') || strcmpi(Method,'FQI')
        [Ave_Payoff,Payoffs,ExBoundary_root] = LSPI_FQI_Performance(Weights,Data,Vol_Data,r,Settings, Ticker, MARKER_plot, method, K);
%         if MARKER_Exer_Boundary
%             ExBoundary_exh = LSPI_FQI_Exercise_Boundary(Weights,Data,mean(mean(Vol_Data)),r,Settings, Ticker, MARKER_plot);
%         end
else    % Method = 'LSM'
       [Ave_Payoff,Payoffs,ExBoundary_root] = LSM_Performance_OriginalPaper(Weights,Data,Vol_Data,r,Settings, Ticker, MARKER_plot, method, K);
       % [Ave_Payoff,Payoffs,ExBoundary_root] = LSM_Performance(Weights,Data,Vol_Data,r,Settings, Ticker, MARKER_plot, method, K);
%         if MARKER_Exer_Boundary
%             ExBoundary_exh = LSM_Exercise_Boundary(Weights,Data,mean(mean(Vol_Data)),r,Settings);
%         end    
end
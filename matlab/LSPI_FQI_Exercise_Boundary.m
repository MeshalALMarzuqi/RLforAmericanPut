function [Boundary] = LSPI_FQI_Exercise_Boundary(Weights,Data,sigma,r,Settings, Ticker, MARKER_plot)            
% Compiling the exercise boundary for the first path. The boundary B_t is solved from the following inequality: 
% max(1 - S_t, 0) = Intrinsic Value > Intrinsic Value Boundary = max(1 - B_t, 0) = Continuation Value(B_t) = Phi_t(B_t)' * Weights
%
Tenor       = size(Data,2) - 1;
Boundary    = zeros(Tenor,1) * NaN;
S_Grid      = 0.1 : 0.001 : 3;
I           = length(S_Grid);
ind         = 1:I;
for t = 1 : Tenor  
    I_Vals          = zeros(I,1) * NaN;
    C_Vals          = zeros(I,1) * NaN;
    for i = 1 : I
        S_t         = S_Grid(i);
        Intr_Value  = max(1 - S_t, 0);         % The intrinsic (exercise) value of a put option with strike 1. 
        Phi_t       = [1; exp(-S_t/2); exp(-S_t/2)*(1-S_t); exp(-S_t/2)*(1-2*S_t+S_t^2/2); sin(-t/Tenor * pi/2 + pi/2); log(Tenor - t); (t/Tenor)^2];
        Cont_Value  = Phi_t' * Weights;  
        I_Vals(i)   = Intr_Value;
        C_Vals(i)   = Intr_Value - Cont_Value; % our f(S) equation
    end  
    flags           = I_Vals > 0 & C_Vals > 0; % (C_Vals > 0) is our K-St-Q>0
    if sum(flags) > 0
        Boundary(t,1) = S_Grid(max(ind(flags)));
    end
end
if MARKER_plot
    figure, plot(Boundary); xlabel('Time[years]'); ylabel('Stock Price [a.u.]'); title(['LSPIorFQI, ' Ticker ', Exercise Boundary, Exhaustive'])
end
end
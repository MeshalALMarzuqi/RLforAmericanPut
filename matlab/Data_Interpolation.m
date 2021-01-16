%% Deals with missing data and interpolates to fill in gaps
function [Mod_Series] = Data_Interpolation(Series)

%
% This function makes sure the series is marked for every time point.
%

T           = length(Series);
MRMI        = zeros(T,1);
NMI         = zeros(T,1);

% Compiling the list of most recent dates when the series was marked.
MRI = 0;
for t = 1:T
   if ~isnan(Series(t))  % check if vector element is not a number and note the time step
        MRI = t;
   end
   MRMI(t) = MRI;
end

% Compiling the list of forthcoming dates when the series will be marked.
NI = 0;
for t = T:-1:2
   if ~isnan(Series(t)) 
        NI = t;
   end
   NMI(t-1) = NI;
end       

% Interpolating between the most recent and forthcoming dates.
Mod_Series = Series;   
for t = 1:T
    if isnan(Mod_Series(t)) && (MRMI(t) > 0)
            Mod_Series(t)    = Mod_Series(MRMI(t));
    end
end
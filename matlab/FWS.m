function [Output] = FWS(Input,Width,Side)

IW      = length(Input);

if Width <= IW
    Output      = Input(1:Width);
else
    Output      = Input;
    
    if strcmpi(Side,'right')
        for i = (IW+1):Width
            Output      = [' ', Output];
        end        
    else
        for i = (IW+1):Width
            Output      = [Output, ' '];
        end
    end
end
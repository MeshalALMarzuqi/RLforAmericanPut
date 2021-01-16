% Bisection Root Searching Algorithm
% e.g. find sqrt of 9
% f(x) = x^2
num = 9; % number to be squarerooted
tol = 0.001;
diff = 99;
count = 0;
low  = 0;
high = num;
while abs(diff) > tol
    mid  = (low + high)/2;
    fom  = mid^2;
    diff = fom - num;
    if diff > 0
        high = mid;
    else
        low = mid;
    end
    count = count + 1;
end
fprintf(['Root of ' num2str(num) ' is ' num2str(mid) '\n'])
fprintf(['Iterated ' num2str(count) ' times.' '\n'])














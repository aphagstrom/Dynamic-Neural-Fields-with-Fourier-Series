function wx = evaluateFourier(x, Wi, Modes)

sum = Wi(1)/2;

for n = 1:Modes-1
    sum = sum + Wi(n+1)*cos(n*x);
end

wx = sum;

 

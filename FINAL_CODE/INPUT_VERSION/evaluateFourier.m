function wx = evaluateFourier(x, Wi_A,Wi_B, Modes)

sum = Wi(1)/2;

for n = 1:Modes-1
    sum = sum + Wi_A(n+1)*cos(n*x)+Wi_B(n+1)*sin(n*x);
end

wx = sum;

 

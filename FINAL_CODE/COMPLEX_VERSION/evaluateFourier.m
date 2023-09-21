function wx = evaluateFourier(x, A, Modes)


sum = 0;

for n = -Modes:Modes
    sum = sum + accessFourierCoeff(n,A,Modes)*exp(1i*n*x);
end


if imag(sum) > 10^(-6)
    crashNow = crashHere
end

wx = real(sum);

return
 

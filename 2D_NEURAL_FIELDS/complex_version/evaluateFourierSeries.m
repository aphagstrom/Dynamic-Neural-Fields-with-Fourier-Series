function fx = evaluateFourierSeries(x,y,ll,h,cmn,Modes)
sum = 0;
for m=-Modes:Modes
   for n = -Modes:Modes
    sum = sum + accessFourierCoeff(m,n,cmn,Modes)*exp(1j*pi*((m*x)./ll+(n*y)./h));
   end
end

if imag(sum) > 10^(-6)
   crashNow = crashHere;
end

fx = real(sum);
return

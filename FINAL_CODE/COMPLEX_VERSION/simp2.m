function T = simp2(f,a,b,N,i)
sum = 0;
numPoints = 10^4;
step = (b-a)/N; 
ex = @(x)(exp(-1i*i*x));

sum=sum+f(1)*ex(a)+f(end)*ex(b)+4*f(5001)*ex(0);
T = step./6*(sum);
return 
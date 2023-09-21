function T = simp1(f,a,b,N,j,i)
sum = 0;
numPoints = 10^4;
step = (b-a)/N; 
ex = @(x)(exp(-1i*i*x));
ex2 = @(x)(exp(1i*j*x));
sum=sum+f(1)*ex(a)*ex2(a)+f(end)*ex(b)*ex2(b)+4*f(5001)*ex(0)*ex2(0);
T = step./6*(sum);
return 
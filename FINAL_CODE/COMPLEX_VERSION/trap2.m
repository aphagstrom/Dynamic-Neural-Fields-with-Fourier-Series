function T = trap2(f,a,b,N,i)
sum = 0;
numPoints = 10^4;
step = (b-a)/N; 
ex = @(x)(exp(-1i*i*x));
x=linspace(-pi,pi,N);
for n = 1:N-2
    sum = sum+2*f(n+1)*ex(x(n+1));
end
sum=sum+f(1)*ex(a)+f(end)*ex(b);
T = step./2*(sum);
return

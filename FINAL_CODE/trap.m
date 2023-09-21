function T = trap(f,a,b,N)

sum = 0;
step = (b-a)/N; 
avg = (f(b)+f(a))/2;
for n = 1:N-1
    xn = a + n*step;
    sum = sum + f(xn);
end

T = step*(avg + sum);
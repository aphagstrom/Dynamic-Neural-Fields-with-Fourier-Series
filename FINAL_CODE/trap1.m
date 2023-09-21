function T = trap1(f,a,b,N)

sum = 0;
step = (b-a)/N; 
avg = (f(1)+f(end))/2;
for n = 1:N-2
    sum = sum + f(n+1);
end

T = step*(avg + sum);
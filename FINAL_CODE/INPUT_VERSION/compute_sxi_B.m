function an = compute_sxi_B(f,n,N,a,b)
%Our integrand is f(x)*cos(nx)
sum = 0;
step = (b-a)/N; 
avg = (f(end)*sin(n*10000)+f(1)*sin(n*1))/2;
for n = 1:N-1
    sum = sum + f(n)*sin(n);
end

T = step*(avg + sum);
an = (1/pi)*(T);
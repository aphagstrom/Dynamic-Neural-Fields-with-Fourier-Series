function an = compute_sxi_A(f,n,N,a,b)
%Our integrand is f(x)*cos(nx)
sum = 0;
step = (b-a)/N; 
avg = (f(end)*cos(n*10000)+f(1)*cos(n*1))/2;
for n = 1:N-1
    sum = sum + f(n)*cos(n);
end

T = step*(avg + sum);
an = (1/pi)*(T);
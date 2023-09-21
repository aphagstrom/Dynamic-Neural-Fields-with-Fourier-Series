function an = compute_wxi_B(f,n,N,a,b)
%Our integrand is f(x)*sin(nx)
sum = 0;
step = (b-a)/N; 
avg = (f(b)*sin(n*b)+f(a)*sin(n*a))/2;
for k = 1:N-1
    xk = a + k*step;
    sum = sum + f(xk)*sin(n*xk);
end

T = step*(avg + sum);
an = (1/pi)*(T);
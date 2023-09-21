function an = compute_wxi_A(f,n,N,a,b)
%Our integrand is f(x)*cos(nx)
sum = 0;
step = (b-a)/N; 
avg = (f(b)*cos(n*b)+f(a)*cos(n*a))/2;
for k = 1:N-1
    xk = a + k*step;
    sum = sum + f(xk)*cos(n*xk);
end

T = step*(avg + sum);
an = (1/pi)*(T);
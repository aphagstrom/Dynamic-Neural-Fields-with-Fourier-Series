function an = compute_complex_wxi(f,n,N,a,b)
% Integrand is: f(x) exp(-inx)
sum = 0;
step = (b-a)/N; 
ex = @(x)(exp(-1i*n*x));
avg = (f(b)*ex(b)+f(a)*ex(a))/2;
for n = 1:N-1
    xn = a + n*step;
    sum = sum + f(xn)*ex(xn);
end
an = step*(avg + sum)/(2 *pi);

return
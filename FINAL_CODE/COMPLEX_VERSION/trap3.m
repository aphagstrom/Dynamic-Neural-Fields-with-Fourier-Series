function T = trap3(f,a,b,N)
sum = 0;
numPoints = 1000;
step = (b-a)/N; 
avg = (f(1)+f(end))/2;
for n = 1:N-1
    the_x1=linspace(-pi,pi,numPoints);
    sum = sum + f(n);
end
T = step*(avg + sum)/(2*pi);
return

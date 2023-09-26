function T2=Trapezoidal_2D(f,N,deltaX,deltaY,x,y,a,b,c,d)


%[X,Y]=meshgrid(x,y);
%meshc(X, Y, f(X,Y))
%grid on
XC_sum=zeros(N,1);
for n = 2:N-1
    XC_sum(n+1) =  XC_sum(n)+f(x(n),c);
end
   XD_sum=zeros(N,1);
for n = 2:N-1
    XD_sum(n+1) = XD_sum(n) + f(x(n),d);
end
    AY_sum= zeros(N,1);
for n = 2:N-1
    AY_sum(n+1) = AY_sum(n) + f(a,y(n));
end
BY_sum= zeros(N,1);
for n = 2:N-1
    BY_sum(n+1) = BY_sum(n) + f(b,y(n));
end
XY_sum= 0;
for i = 2:N
    for j=2:N
    XY_sum = XY_sum + f(x(j),y(i));
    end
end
T2=deltaX.*deltaY./4*(f(a,c)+f(a,d)+f(b,c)+f(b,d)+2*XC_sum(end)+2*XD_sum(end)+2*AY_sum(end)+2*BY_sum(end)+4*XY_sum(end));
end
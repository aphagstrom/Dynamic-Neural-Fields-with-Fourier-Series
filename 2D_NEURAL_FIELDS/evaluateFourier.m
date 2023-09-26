function fx = evaluateFourier(x,y,amn,bmn,cmn,dmn, Modes,l,h,lambdamn)
sum = 0;
for m = 0:Modes-1
    for n=0:Modes-1
       sum = sum + lambdamn(m+1,n+1).*(amn(m+1,n+1).*cos(pi*m*x./l).*cos(pi*n*y./h)+bmn(m+1,n+1).*sin(pi*m*x./l).*cos(pi*n*y./h)+cmn(m+1,n+1).*cos(pi*m*x./l).*sin(pi*n*y./h)+dmn(m+1,n+1).*sin(pi*m*x./l).*sin(pi*n*y./h));
    end
end
fx = sum;

 

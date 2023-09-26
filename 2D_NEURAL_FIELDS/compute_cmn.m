function [cmn]=compute_cmn(W,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h)
XC_sum=0;
for i=1:N-1
  for j=1:N-1
    XC_sum=  XC_sum+W(x(i),c)*cos(pi*m*x(i)./l)*sin(pi*n*y(j)./h);
  end
end
XD_sum=0;
for i=1:N-1
   for j=1:N-1
    XD_sum = XD_sum + W(x(i),d)*cos(pi*m*x(i)./l)*sin(pi*n*y(j)./h);
   end
end
AY_sum= 0;
for i=1:N-1
  for j=1:N-1
    AY_sum = AY_sum + W(a,y(j))*cos(pi*m*x(i)./l)*sin(pi*n*y(j)./h);
  end
end
BY_sum= 0;
for i=1:N-1
  for j=1:N-1
      BY_sum = BY_sum + W(b,y(j))*cos(pi*m*x(i)./l)*sin(pi*n*y(j)./h);
  end
end
XY_sum=0;
for i = 1:N-1
  for j=1:N-1
     XY_sum = XY_sum + W(x(i),y(j))*cos(pi*m*x(i)./l)*sin(pi*n*y(j)./h);
   end
end
cmn=deltaX.*deltaY./4*(W(a,c)+W(a,d)+W(b,c)+W(b,d)+2*XC_sum+2*XD_sum+2*AY_sum+2*BY_sum+4*XY_sum);
end

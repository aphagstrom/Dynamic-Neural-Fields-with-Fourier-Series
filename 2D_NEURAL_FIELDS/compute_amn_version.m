function [amn1]=compute_amn_version(W,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h)
XC_sum=0;
for i=1:N-1
   for j=1:N-1
      XC_sum=XC_sum+W(i,1)*cos(pi*m*x(i)./l)*cos(pi*n*y(j)./h);
   end
end
XD_sum=0;
for i=1:N-1
   for j=1:N-1
    XD_sum = XD_sum + W(i,end)*cos(pi*m*x(i)./l)*cos(pi*n*y(j)./h);
   end
end
AY_sum= 0;
for i=1:N-1
  for j=1:N-1
    AY_sum = AY_sum + W(1,j)*cos(pi*m*x(i)./l)*cos(pi*n*y(j)./h);
  end
end
BY_sum=0;
for i=1:N-1
  for j=1:N-1
      BY_sum = BY_sum + W(end,j)*cos(pi*m*x(i)./l)*cos(pi*n*y(j)./h);
  end
end
XY_sum=0;
for i = 1:N-1
    for j=1:N-1
      XY_sum = XY_sum + W(i,j)*cos(pi*m*x(i)./l)*cos(pi*n*y(j)./h);
    end
end
amn1=deltaX.*deltaY./4*(W(1,1)+W(1,end)+W(end,1)+W(end,end)+2*XC_sum+2*XD_sum+2*AY_sum+2*BY_sum+4*XY_sum);
end

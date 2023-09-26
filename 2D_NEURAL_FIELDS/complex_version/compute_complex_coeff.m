function [cn]=compute_complex_coeff(W,N,deltaX,deltaY,x,y,m,n,ll,h,a,b,c,d)

ex = @(x,y)(exp((-1j*pi)*((m*x)./ll+(n*y)./h)));

XC_sum=0;
for i=2:N-1
      XC_sum=XC_sum+W(x(i),c)*ex(x(i),c); 
end
XD_sum=0;
for i=2:N-1
    XD_sum = XD_sum + W(x(i),d)*ex(x(i),d);
end
AY_sum= 0;
  for j=2:N-1
    AY_sum = AY_sum + W(a,y(j))*ex(a,y(j));
  end

BY_sum=0;
  for j=2:N-1
      BY_sum = BY_sum + W(b,y(j))*ex(b,y(j));
  end
XY_sum=0;
for  i= 2:N-1
    for j=2:N-1
      XY_sum = XY_sum + W(x(i),y(j))*ex(x(i),y(j));
    end
end
cn=deltaX.*deltaY./4*(W(a,c)*ex(a,c)+W(a,d)*ex(a,d)+W(b,c)*ex(b,c)+W(b,d)*ex(b,d)+2*XC_sum+2*XD_sum+2*AY_sum+2*BY_sum+4*XY_sum);
end

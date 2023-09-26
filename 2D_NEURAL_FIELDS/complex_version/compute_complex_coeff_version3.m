function [cn]=compute_complex_coeff_version3(W,N,deltaX,deltaY,m,n,ll,h,x,y,a,b,c,d,aa,bb)
ex = @(x,y)(exp((-1j*pi)*((m*x)./ll+(n*y)./h)));
ex2 = @(x,y)(exp((1j*pi)*((aa*x)./ll+(bb*y)./h)));

XC_sum=0;
for i=2:N-1

      XC_sum=XC_sum+W(i,1).*ex(x(i),c).*ex2(x(i),c);

end
XD_sum=0;
for i=2:N-1
    XD_sum = XD_sum + W(i,end).*ex(x(i),d).*ex2(x(i),d);

end
AY_sum= 0;
  for j=2:N-1
    AY_sum = AY_sum + W(1,j).*ex(a,y(j)).*ex2(a,y(j));
  end

BY_sum=0;

  for j=2:N-1
      BY_sum = BY_sum + W(end,j).*ex(b,y(j)).*ex2(b,y(j));
  end

XY_sum=0;
for i = 2:N-1
    for j=2:N-1
      XY_sum = XY_sum + W(i,j).*ex(x(i),y(j)).*ex2(x(i),y(j));
    end
end
cn=deltaX.*deltaY./4*(W(1,1)*ex(a,c)*ex2(a,c)+W(1,end)*ex(a,d)*ex2(a,d)+W(end,1)*ex(b,c)*ex2(b,c)+W(end,end)*ex(b,d)*ex2(b,d)+2*XC_sum+2*XD_sum+2*AY_sum+2*BY_sum+4*XY_sum);
end

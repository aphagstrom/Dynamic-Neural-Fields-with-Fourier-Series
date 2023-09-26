function [cn]=compute_complex_coeff_version2(W,N,deltaX,deltaY,m,n,ll,h,x,y,a,b,c,d)

ex = @(x,y)(exp((-1i*pi)*((m*x)./ll+(n*y)./h)));

XC_sum=0;
for i=1:N-1
   %for j=1:N-1
      XC_sum=XC_sum+W(i,1,:,:)*ex(x(i),c);
   %end
end
XD_sum=0;
for i=1:N-1
  % for j=1:N-1
    XD_sum = XD_sum + W(i,end,:,:)*ex(x(i),d);
   %end
end
AY_sum= 0;
%for i=1:N-1
  for j=1:N-1
    AY_sum= AY_sum + W(1,j,:,:)*ex(a,y(j));
  end
%end
BY_sum=0;
%for i=1:N-1
  for j=1:N-1
      BY_sum = BY_sum + W(end,j,:,:)*ex(b,y(j));
  end
%end
XY_sum=0;
for i = 1:N-1
    for j=1:N-1
      XY_sum = XY_sum + W(i,j,:,:)*ex(x(i),y(j));
    end
end
cn=deltaX.*deltaY./4*(W(1,1,:,:)*ex(a,a)+W(1,end,:,:)*ex(a,d)+W(end,1,:,:)*ex(b,c)+W(end,end,:,:)*ex(b,d)+2*XC_sum+2*XD_sum+2*AY_sum+2*BY_sum+4*XY_sum);
end

function S=simps4(f,N,m,n,ll,h,aa,bb)




ex = @(x,y)(exp((-1j*pi)*((m*x)./ll+(n*y)./h)));
ex2 = @(x,y)(exp((1j*pi)*((aa*x)./ll+(bb*y)./h)));





S=4*pi^2./36*(16*f(N./2+1,N./2+1)*ex(0,0)*ex2(0,0)+4*f(N./2+1,end)*ex(0,pi)*ex2(0,pi)+4*f(N./2+1,1)*ex(0,-pi)*ex2(0,-pi)+4*f(end,N./2+1)*ex(pi,0)*ex2(pi,0)+...
    f(end,end)*ex(pi,pi)*ex2(pi,pi)+f(end,1)*ex(pi,-pi)*ex2(pi,-pi)+4*f(1,N./2+1)*ex(-pi,0)*ex2(-pi,0)+f(1,end)*ex(-pi,pi)*ex2(-pi,pi)+f(1,1)*ex(-pi,-pi)*ex2(-pi,-pi));


end

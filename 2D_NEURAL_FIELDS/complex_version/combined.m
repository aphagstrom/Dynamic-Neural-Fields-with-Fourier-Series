function cmn_new  = setFourierCoef(m,n,cmn,cmn_new,Modes)
cmn_new(Modes+1+m,Modes+1+n) = cmn;
return
end

function [F] = neuralField2D_F_complex(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d)
Fi = zeros(2*Modes+1, 2*Modes+1);
S=@(x) 1./(1+exp(-Beta.*((x-hh))));
for m = -Modes:Modes
    for n=-Modes:Modes
    for k=1:N
        for l=1:N
                
                   thisEval(k,l) = evaluateFourierSeries(x(k),y(l),ll,h,dmn,Modes);
                   thisSEval(k,l) = S(thisEval(k,l));
        end  
    end
                   Integral_1(Modes+m+1,Modes+n+1)=compute_complex_coeff_version(thisSEval,N,deltaX,deltaY,m,n,ll,h,x,y,a,b,c,d);
                   Fi(Modes+m+1,Modes+n+1)=cmn(Modes+m+1,Modes+n+1)*Integral_1(Modes+m+1,Modes+n+1)-dmn(Modes+m+1,Modes+n+1);
     end
end
F = Fi;
return
end

function [DF,thisSEval] = neuralField2D_DF(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d)

DFi = zeros(2*Modes+1, 2*Modes+1,2*Modes+1,2*Modes+1);               %DFi

Sprime=@(x) Beta.*exp(-Beta.*(x-hh))./(1+exp(-Beta.*((x-hh)))).^2;    %S'(x)

for m = -Modes:Modes
    for n=-Modes:Modes
       for aa=-Modes:Modes
           for bb=-Modes:Modes
                    for k=1:N
                        for l=1:N
              
                                 thisEval(k,l) = evaluateFourierSeries(x(k),y(l),ll,h,dmn,Modes);
                                 thisSEval(k,l) = Sprime(thisEval(k,l));
                                 %thisSEval(k,l,Modes+aa+1,Modes+bb+1) = Sprime(thisEval(k,l))*exp((1i*pi)*((aa*x(k))./ll+(bb*y(l))./h));
                          
                                 %thisSEval2(k,l,Modes+1+aa,Modes+1+bb)=thisSEval.*exp((1i*pi)*((aa*x(k))./ll+(bb*y(l))./h));
                                
                        end
                    
                    end
                              
                                
                                %thisSEval2=squeeze(thisSEval2(:,:,Modes+1+aa,Modes+1+bb));
                               
                                 
                                 if (m==aa) && (n==bb)

                                 DFi(Modes+m+1,Modes+n+1,Modes+aa+1,Modes+bb+1)=cmn(Modes+m+1,Modes+n+1)*compute_complex_coeff_version3(thisSEval,N,deltaX,deltaY,m,n,ll,h,x,y,a,b,c,d,aa,bb)-1;       
       
                                 else 

                                 DFi(Modes+m+1,Modes+n+1,Modes+aa+1,Modes+bb+1)= cmn(Modes+m+1,Modes+n+1)*compute_complex_coeff_version3(thisSEval,N,deltaX,deltaY,m,n,ll,h,x,y,a,b,c,d,aa,bb);
                                 
                                 end
           end
       end
            
    end
end

DF=DFi;
return
end


function A = emptyFourier(Modes)
A = zeros(2*Modes+1, 2*Modes+1);
return
end

function fx = evaluateFourierSeries(x,y,ll,h,cmn,Modes)
sum = 0;
for m=-Modes:Modes
   for n = -Modes:Modes
    sum = sum + accessFourierCoeff(m,n,cmn,Modes)*exp(1i*pi*((m*x)./ll+(n*y)./h));
   end
end

if imag(sum) > 10^(-6)
   crashNow = crashHere;
end
fx = real(sum);
return
end

function cn = computeComplexFourierCoefficients(f)
% Integrand is: f(x) exp(-inx)
N = 10^5;
a = -pi;
b = pi;
sum = 0;
step = (b-a)/N; 

ex = @(x)(exp(-1i*n*x));

avg = (f(b)*ex(b)+f(a)*ex(a))/2;
for n = 1:N-1
    xn = a + n*step;
    sum = sum + f(xn)*ex(xn);
end
cn = step*(avg + sum)/(2*pi);
return
end

function [cn]=compute_complex_coeff_version3(W,N,deltaX,deltaY,m,n,l,h,x,y,a,b,c,d,aa,bb)

ex = @(x,y)(exp((-1i*pi)*((m*x)./l+(n*y)./h)));
ex2 = @(x,y)(exp((1i*pi)*((aa*x)./l+(bb*y)./h)));


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

function [cn]=compute_complex_coeff_version2(W,N,deltaX,deltaY,m,n,l,h,x,y,a,b,c,d)

ex = @(x,y)(exp((-1i*pi)*((m*x)./l+(n*y)./h)));

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
return
end

function [cn]=compute_complex_coeff_version(W,N,deltaX,deltaY,m,n,l,h,x,y,a,b,c,d)

ex = @(x,y)(exp((-1i*pi)*((m*x)./l+(n*y)./h)));

XC_sum=0;
for i=2:N-1
 
      XC_sum=XC_sum+W(i,1)*ex(x(i),c);

end
XD_sum=0;
for i=2:N-1

    XD_sum = XD_sum + W(i,end)*ex(x(i),d);
 
end
AY_sum= 0;

 for j=2:N-1
    AY_sum = AY_sum + W(1,j)*ex(a,y(j));
 end

BY_sum=0;

  for j=2:N-1
      BY_sum = BY_sum + W(end,j)*ex(b,y(j));
  end
XY_sum=0;
for i = 2:N-1
    for j=2:N-1
      XY_sum = XY_sum + W(i,j)*ex(x(i),y(j));
    end
end
cn=deltaX.*deltaY./4*(W(1,1)*ex(a,c)+W(1,end)*ex(a,d)+W(end,1)*ex(b,c)+W(end,end)*ex(b,d)+2*XC_sum+2*XD_sum+2*AY_sum+2*BY_sum+4*XY_sum);
return
end
function a = accessFourierCoeff(m,n,cmn,Modes)
if (-Modes <= n) && (n <= Modes) && (-Modes <= m) && (m <= Modes)
    a = cmn(Modes+1+m,Modes+1+n);
else
    a = zeros(size(cmn));
end
return
end

function [cn]=compute_complex_coeff(W,N,deltaX,deltaY,x,y,m,n,ll,h,a,b,c,d)
ex = @(x,y)(exp((-1i*pi)*((m*x)./ll+(n*y)./h)));

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
return
end











W=@(x,y) exp(-(x.^2+y.^2)-0.17*(exp(-0.2*(x.^2+y.^2)))) ; %function handle for 2D connectivity 

 %%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
Beta=10;
hh=0.48;
ll=pi;
h=pi;
a=-ll; %x left endpoint
b=ll;  %x right endpoint
c=-h; %y left endpoint
d=h;  %y right endpoint
%N=300; %num ticks in x and y directions
%x=linspace(-pi,pi,N);
%y=linspace(-pi,pi,N);
%deltaX=length(x)./N;
%deltaY=length(y)./N;
%%%%%%%%%%%%%%%%% GRID %%%%%%%%%%%%%%%%%%%%%%%
deltaX=(ll-(-ll))./60;  %grid spacing in x-direction
deltaY=(h-(-h))./60;  %grid spacing in y-direction
x = a:deltaX:b;   %x grid
y = c:deltaY:d;   % y grid
N=(b-a)./deltaX;  %number of points in grid
Modes=1;     %number of modes
%%%%%%%%%%%%%%%%%%%%%%%% Complex Version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmn=emptyFourier(Modes);
for m=-Modes:Modes
   for n=-Modes:Modes
       n
       cmn(Modes+1+m,Modes+1+n)=1./(2*ll*2*h)*compute_complex_coeff(W,N,deltaX,deltaY,x,y,m,n,ll,h,a,b,c,d);
   end
end

W_values=zeros(N,N);
for k=1:N
    for j=1:N    
        Wx(k,j)=evaluateFourierSeries(x(k),y(j),ll,h,cmn,Modes);
        W_values(k,j)=W(x(k),y(j));
    end
end

dmn=cmn./2; %Ui coefficients
figure(1)
surf(Wx);
figure(2)
surf(W_values);


%% Error
D = abs(Wx-W_values).^2;
MSE_W = sum(D(:))/numel(Wx); %MSE Error of W
MSE = norm((W_values-Wx),'fro')^2/numel(Wx); %Frobenius norm


%% Calculate F and DF
F = neuralField2D_F_complex(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d);
DF = neuralField2D_DF(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d);
DF=reshape(DF,(2*Modes+1)^2,(2*Modes+1)^2,1);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVE THE SYSTEM OF ODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DFpartial = neural2D_partial(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,m,n,aa,bb);


tau=10^(-7);
delta0_F=dmn;
delta0_F(3,2)=dmn(3,2)+tau;
delta0_F=neuralField2D_F_complex(cmn,delta0_F,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d);
partial1=(delta0_F(3,2)-F(3,2))./tau;




return
tau=10^(-7);
delta0_F=dmn;
delta0_F(5,1)=dmn(5,1)+tau;
delta0_F=neural2D_F(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,1,1,a,b,c,d);
partial2=(delta0_F-F(5,1))./tau;

%NEWTON'S METHOD
X=zeros(size(F));
X=reshape(X,(2*Modes+1)^2,1);
U=zeros(N,N);
hold on
F=reshape(F,(2*Modes+1)^2,1);
DF=reshape(DF,(2*Modes+1)^2,(2*Modes+1)^2);

%FINITE DIFFERENCING TEST
DF_test = zeros((2*Modes+1)^2,(2*Modes+1)^2);
for index = 1:(2*Modes+1)^2
    tau = 10^(-7);
    delta_Ui=reshape(dmn,[(2*Modes+1)^2,1]);
    delta_Ui(index) = delta_Ui(index) + tau;
    delta_Ui=reshape(delta_Ui,[2*Modes+1,2*Modes+1]);
    delta0_F = neuralField2D_F_complex(cmn,delta_Ui,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d);
    delta0_F=reshape(delta0_F,(2*Modes+1)^2,1);
    partial1 = (delta0_F - F)/tau;
    DF_test(:,index) = partial1;
end


%FRECHET DERIVATIVE TEST
 H = zeros((2*Modes+1)^2, 1);
 scale = 10^(-7);
 for k = 1:(2*Modes+1)^2
     H(k) = scale*(0.25^(k-1))*rand(1,1);
 end
 dmn=reshape(dmn,(2*Modes+1)^2,1)+H;
 dmn=reshape(dmn,(2*Modes+1),(2*Modes+1));
 Fdelta = neuralField2D_F_complex(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d);
 Fdelta=reshape(Fdelta,(2*Modes+1)^2,1);
 error_Frechet = norm(Fdelta - F - DF*H, "inf");



%check 10 columns
colerror_1=norm(DF_test(:,1)-DF(:,1));
colerror_2=norm(DF_test(:,2)-DF(:,2));
colerror_3=norm(DF_test(:,3)-DF(:,3));
colerror_4=norm(DF_test(:,4)-DF(:,4));
colerror_5=norm(DF_test(:,5)-DF(:,5));
colerror_6=norm(DF_test(:,6)-DF(:,6));
colerror_7=norm(DF_test(:,7)-DF(:,7));
colerror_8=norm(DF_test(:,8)-DF(:,8));
colerror_9=norm(DF_test(:,9)-DF(:,9));

return



YY=DF\F;

hold on

for t=1:20
                    X=X-YY; 
                    X=reshape(X,[2*Modes+1,2*Modes+1]);
                    %%%%%%%%%%%% Calculate the new Fourier coefficients %%%%%%%%%%%%%%%
                    F = neuralField2D_F_complex(cmn,X,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d);
                   
                    DF= neuralField2D_DF(cmn,X,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d);
                    F=reshape(F,(2*Modes+1)^2,1);
                     error=norm(F,inf)
                    DF=reshape(DF,(2*Modes+1)^2,(2*Modes+1)^2,1,1);
                    DF=squeeze(DF);
                    YY=DF\F;
                    X=reshape(X,[(2*Modes+1)^2,1]);
                    X=X-YY; 
                    X=reshape(X,[2*Modes+1,2*Modes+1]);
                    
                   for k=1:N
                       for l=1:N
                        U(k,l)=evaluateFourierSeries(x(k),y(l),ll,h,X,Modes);
                       end
                   end
                   X=reshape(X,[(2*Modes+1)^2,1]);
                   figure(3)
                   imagesc(U);
                   
                 
        
end

error_Newton=norm(F,inf) 
   %DF(:,:,2:end,:) = [];
                  % DF(:,:,:,2:end) = [];
                   %DF= squeeze(DF); 







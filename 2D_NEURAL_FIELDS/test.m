W=@(x,y) exp(-(x.^2+y.^2)-0.17*(exp(-0.2*(x.^2+y.^2)))) ; %function handle for 2D connectivity 
l=pi;
h=pi;
a=-l; %x left endpoint
b=l;  %x right endpoint
c=-h; %y left endpoint
d=h;  %y right endpoint
deltaX=(l-(-l))./60;  %grid spacing in x-direction
deltaY=(h-(-h))./60;  %grid spacing in y-direction
x = a:deltaX:b;   %x grid
y = c:deltaY:d;   % y grid
N=(b-a)./deltaX;  %number of points in grid
Modes=10;     %number of modes

%T2=Trapezoidal_2D(f,N,deltaX,deltaY,x,y,a,b,c,d);
amn=zeros(Modes,Modes); % amn coefficient
for m=0:Modes-1
    for n=0:Modes-1
         amn(m+1,n+1)=1./(l*h)*compute_amn(W,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h); %Compute Fourier Coefficients amn
    end
end
bmn=zeros(Modes,Modes); %bmn coefficient
for m=0:Modes-1
    for n=0:Modes-1
         bmn(m+1,n+1)=1./(l*h)*compute_bmn(W,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h); %Compute Fourier Coefficients amn
    end
end
cmn=zeros(Modes,Modes); %cmn coefficient
for m=0:Modes-1
    for n=0:Modes-1
         cmn(m+1,n+1)=1./(l*h)*compute_cmn(W,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h); %Compute Fourier Coefficients amn
    end
end
dmn=zeros(Modes,Modes);  %dmn coefficient
for m=0:Modes-1
    for n=0:Modes-1
       dmn(m+1,n+1)=1./(l*h)*compute_dmn(W,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h); %Compute Fourier Coefficients amn
    end
end


for m=0:Modes-1      %lambda function used in evaluating Fourier Series
    for n=0:Modes-1    
       if (m==0 && n==0)
    lambdamn(m+1,n+1)=1./4;
       elseif (m>0 && n==0)
    lambdamn(m+1,n+1)=1./2;
       elseif (m==0 && n>0)
    lambdamn(m+1,n+1)=1./2;
       elseif (m>0 && n>0) 
    lambdamn(m+1,n+1)=1;
       end
    end
end


Wx=zeros(N,N);
W_values=zeros(N,N);


for i=1:N
    for j=1:N
       Wx(i,j) = evaluateFourier(x(i),y(j),amn,bmn,cmn,dmn,Modes,l,h,lambdamn); %W_values from Fourier series evaluation
       W_values(i,j)=W(x(i),y(j));
    end
end

%MSE Error of W
D = abs(Wx-W_values).^2;
MSE_W = sum(D(:))/numel(Wx);

%Plot of W derived from Fourier Series
figure(1)
imagesc(Wx);

%Plot of W derived from function handle
figure(2)
imagesc(W_values);



%%




%%%%%%%%%%%%%%%%%%%%%% COMPUTING F %%%%%%%%%%%%%%%%%%%%%%

Beta=20;
h=1./3;   
emn=amn./3;
fmn=bmn./3;
gmn=cmn./3;
hmn=dmn./3;

[F1i,F2i,F3i,F4i] = neuralField2D_F(amn,bmn,cmn,dmn,emn,fmn,gmn,hmn,lambdamn,Beta,h,N,deltaX,deltaY,x,y,a,b,c,d,Modes,l);


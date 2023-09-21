clear vars
format long
%Input
rgbImage = imread("cameraman.tif");
I = rgbImage(1:100,1:100);
I=reshape(I,[10000,1]);
%I=I./10000;
I=rescale(I, 0, 1);

%%%%%%%%%%%%%%%  COMPUTE W(X) WITH FOURIER COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Modes=30;
N=10^4;
%functions
sigma=0.0350195;
A=1.2;
B=0.2;
gamma=2;
g=@(x) exp(-abs(x).^2./sigma.^2);
ggamma=@(x) exp(-abs(x).^2./(gamma*sigma).^2);
%W=@(x)(10.*exp(-4.*x.^2)-6.*exp(-x.^2)); %W(x);
W=@(x) A*g(x)-B*ggamma(x); %W(x)
a1=-pi;   %Lower Bound
b1=pi;    %Upper Bound
%w0=(1./(pi)).*trap(W,a1,b1,N); %W_0
Wi_A=zeros(1,Modes); 
Wi_B=zeros(1,Modes);%initialize Wi Fourier Coefficients
Si_A=zeros(1,Modes);
Si_B=zeros(1,Modes);
for n=0:Modes-1
Wi_A(n+1)=compute_wxi_A(W,n,N,a1,b1); %Compute Fourier Coefficients for W(x)
Wi_B(n+1)=compute_wxi_B(W,n,N,a1,b1); %Compute Fourier Coefficients for W(x)
Si_A(n+1)=compute_sxi_A(I,n,N,a1,b1); %Compute Fourier Coefficients for W(x)
Si_B(n+1)=compute_sxi_B(I,n,N,a1,b1);
end
numPoints=10000;
the_x1=linspace(-pi,pi,numPoints);
the_values1=zeros(1,numPoints);
the_fourierSeries1=zeros(1,numPoints);
the_fourierSeriesS=zeros(1,numPoints);
for k=1:numPoints
    this_x1=the_x1(k);
    the_values1(k)=W(this_x1);
    %the_fourierSeries1(k)=evaluateFourier(this_x1,Wi,Modes);
    the_fourierSeriesS(k)=evaluateFourier(this_x1,Si_A,Si_B,odes);
end
figure
plot(the_x1,the_fourierSeries1, 'o');
hold on
plot(the_x1,the_values1,'o');
hold on
plot(the_x1,the_fourierSeriesS,'o');

x = 2;
error_W = abs(W(x)-evaluateFourier(x,Wi,Modes));
%%
%%%%%%%%%%%%%%%%%%%%%%% COMPUTE F AND DF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1;
b=2;
Ui_A=0*Wi;
Ui_B=0*Wi;
%call Fi
U0=zeros(1,Modes);
U = U0;
dt = 0.02;
%options = odeset('OutputFcn',@odetpbar,'RelTol',1e-3,'AbsTol',1e-3);
%tspan = (0:dt:2);
U0 = U(:,end);
S=@(x) 1./(1+exp(-a.*x+b));
eta=0.001;
mu=0.1;

F=I_neuralField_F_w_input(Ui, Wi, Modes, a,b,Si)

%for t=1:5
%[tout,Unew] = ode45(@(t,U) I_neuralField_F_w_input(Ui, Wi, Modes, a,b,Si),tspan,-Ui,options);
%delta_b=dt*(eta*(1-2+1./mu)*S(max(Unew(:,1)))+1./mu*S(U_new(:,1)).^2);
%delta_a=dt*(eta./a+Unew(max(S)*delta_b);
%a=a+delta_a;
%b=b+delta_b;
%end
Unew=Ui';

%EULER's METHOD
for i=1:2/dt
Unew(:,i+1)=Unew(:,i)+dt*I_neuralField_F_w_input(Unew(:,i), Wi_A, Wi_B, Modes, a,b,Si);
%Ui=Unew;

%delta_b=dt*(eta*(1-2+1./mu).*max(S((Unew(:,i))),[],1)+1./mu.*max(S((Unew(:,i))),[],1).^2);
%delta_a=dt*(eta./a+S(max(S((Unew(:,i))),[],1))*delta_b);
%a=a+delta_a;
%b=b+delta_b;




end


%%
for tt=1:2/dt
    for k=1:numPoints
        this_x1=the_x1(k);
        UU(k,tt)=evaluateFourier(this_x1,Unew(:,tt),Modes);
        SS(k)=evaluateFourier(this_x1,Si_A,Si_B,Modes);
    end
end
figure
hold on
%%
%for tt=1:10./dt
plot(the_x1,UU(:,end),'r')
plot(the_x1,UU(:,end-5),'b')
plot(the_x1,SS,'g')
%end


return

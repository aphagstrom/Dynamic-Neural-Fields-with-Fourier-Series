clear vars
format long
%%%%%%%%%%%%%%%  COMPUTE W(X) WITH FOURIER COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Modes=5;
N=10^4;
%functions
%a = 0;
%b = 3;
g_ex=0.1;
g_inh=0.3;
sigma_ex=.1;
sigma_inh=0.2;
%W=@(x) g_ex./(sqrt(2.*pi).*sigma_ex).*exp(-x^2./(2*sigma_ex).^2)-g_inh./(sqrt(2*pi).*sigma_inh).*exp(-x.^2./(2*sigma_inh.^2)); %W(x)
W=@(x)(10.*exp(-4.*x.^2)-6.*exp(-x.^2)); 



a1=-pi;  %Lower Bound
b1=pi;   %Upper Bound



%w0=(1./(pi)).*trap(W,a1,b1,N); %W_0
Wi=zeros(2*Modes+1,1);               %initialize Wi Fourier Coefficients
for n=-Modes:Modes
Wi(Modes+1+n)=compute_complex_wxi(W,n,N,a1,b1); %Compute Fourier Coefficients for W(x)
end
numPoints=10.^4;
the_x1=linspace(-pi,pi,numPoints);
the_values1=zeros(numPoints,1);
the_fourierSeries1=zeros(numPoints,1);
for k=1:numPoints
    this_x1=the_x1(k);
    the_values1(k)=W(this_x1);
    the_fourierSeries1(k)=evaluateFourier(this_x1,Wi,Modes);
end
figure
plot(the_x1,the_fourierSeries1, 'o');
hold on
plot(the_x1,the_values1,'o');

%%
%%%%%%%%%%%%%%%%%%%%%%% COMPUTE F AND DF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Beta=20;
h=1./3;   
Ui=3*Wi;
dt=0.01;
%call Fi
F = neuralField_F(Ui, Wi, Modes, Beta, h);
%call DFi
DF = neuralfield_DF(Ui, Wi,Modes,Beta,h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINITE DIFFERENCING TEST
DF_test = zeros(2*Modes+1, 2*Modes+1);
for index = -Modes:Modes
    tau = 10^(-7);
    delta_Ui = Ui;
    delta_Ui(index+Modes+1) = Ui(index+Modes+1) + tau;
    delta0_F =  neuralField_F(delta_Ui, Wi, Modes, Beta, h);
    partial1 = (delta0_F - F)/tau;
    DF_test(:, index+Modes+1) = partial1;
end


%check 10 columns
colerror_1=norm(DF_test(:,1)-DF(:,1));
colerror_2=norm(DF_test(:,2)-DF(:,2));
colerror_3=norm(DF_test(:,3)-DF(:,3));
colerror_4=norm(DF_test(:,4)-DF(:,4));
colerror_5=norm(DF_test(:,5)-DF(:,5));
%colerror_6=norm(DF_test(:,6)-DF(:,6));
%colerror_7=norm(DF_test(:,7)-DF(:,7));
%colerror_8=norm(DF_test(:,8)-DF(:,8));
%colerror_9=norm(DF_test(:,9)-DF(:,9));
%colerror_10=norm(DF_test(:,10)-DF(:,10));
error_finitedifferencing=norm(DF-DF_test,inf);


%Unew=3*Ui;
%options = odeset('OutputFcn',@odetpbar,'RelTol',1e-3,'AbsTol',1e-3);
%tspan = (0:dt:2);
%[tout,Unew] = ode45(@(t,U) I_neuralField_F_w_input(Ui, Wi, Modes, a,b,Si),tspan,-Ui,options);




%Unew=Ui';
Unew=3*Ui;
%EULER's METHOD
for i=1:10/dt
Unew(:,i+1)=Unew(:,i)+dt*neuralField_F(Unew(:,i), Wi, Modes, Beta,h);
end
%%
for tt=1:10/dt
    for k=1:numPoints
        this_x1=the_x1(k);
        UU(k,tt)=evaluateFourier(this_x1,Unew(:,tt),Modes);
    end
end
figure
hold on
%%
for tt=1:10./dt
plot(the_x1,UU(:,tt),'k')
pause(.1)
end



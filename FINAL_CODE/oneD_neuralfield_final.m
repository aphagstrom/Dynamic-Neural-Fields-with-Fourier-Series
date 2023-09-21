clear vars
format long
%%%%%%%%%%%%%%%  COMPUTE W(X) WITH FOURIER COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Modes=10;
N=10^4;
%functions
W=@(x)(10.*exp(-4.*x.^2)-6.*exp(-x.^2)); %W(x)
a1=-pi;   %Lower Bound
b1=pi;    %Upper Bound
%w0=(1./(pi)).*trap(W,a1,b1,N); %W_0
Wi=zeros(1,Modes);               %initialize Wi Fourier Coefficients
for n=0:Modes-1
Wi(n+1)=compute_wxi(W,n,N,a1,b1); %Compute Fourier Coefficients for W(x)
end
numPoints=100;
the_x1=linspace(-pi,pi,numPoints);
the_values1=zeros(1,numPoints);
the_fourierSeries1=zeros(1,numPoints);
for k=1:numPoints
    this_x1=the_x1(k);
    the_values1(k)=W(this_x1);
    the_fourierSeries1(k)=evaluateFourier(this_x1,Wi,Modes);
end
figure
plot(the_x1,the_fourierSeries1, 'o');
hold on
plot(the_x1,the_values1,'o');

x = 2;
error_W = abs(W(x)-evaluateFourier(x,Wi,Modes));
%%
%%%%%%%%%%%%%%%%%%%%%%% COMPUTE F AND DF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Beta=20;
h=1./3;   
Ui=2*Wi;
%call Fi
F = neuralField_F(Ui, Wi, Modes, Beta, h);

%call DFi
DF = neuralfield_DF(Ui, Wi,Modes,Beta,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINITE DIFFERENCING TEST
DF_test = zeros(Modes, Modes);
for index = 1:Modes
    tau = 10^(-7);
    delta_Ui = Ui;
    delta_Ui(index) = Ui(index) + tau;
    delta0_F =  neuralField_F(delta_Ui, Wi, Modes, Beta, h);
    partial1 = (delta0_F - F)/tau;
    DF_test(:, index) = partial1;
end

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
colerror_10=norm(DF_test(:,10)-DF(:,10));

error_finitedifferencing=norm(DF-DF_test,inf);
%FRECHET DERIVATIVE TEST
 H = zeros(Modes, 1);
 scale = 10^(-7);
 for k = 1:Modes
     H(k) = scale*(0.25^(k-1))*rand(1,1);
 end
 Fdelta = neuralField_F(Ui+H', Wi, Modes, Beta, h);
error_Frechet = norm(Fdelta - F - DF*H, "inf");


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVE THE SYSTEM OF ODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NEWTON'S METHOD
YY=-DF\F;
X=Ui';
for t=1:10  %Why do we need a "t" here? and why does YY eventually become 0
    X=X+YY;    
    F = neuralField_F(X, Wi, Modes, Beta, h);
    error=norm(F,inf)
    DF= neuralfield_DF(X, Wi,Modes,Beta,h);
    YY=-DF\F;
end
X=X+YY;    
F = neuralField_F(X, Wi, Modes, Beta, h);
error_Newton=norm(F,inf)   %What exactly does this do.
the_x2=linspace(-pi,pi,numPoints);
the_fourierSeries2=zeros(numPoints,1);
for k=1:numPoints
    this_x2=the_x1(k);
    the_fourierSeries2(k)=evaluateFourier(this_x2,X,Modes);
end
plot(the_x2,the_fourierSeries2,'o');
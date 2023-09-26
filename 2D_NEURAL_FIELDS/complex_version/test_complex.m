W=@(x,y) exp(-(x.^2+y.^2)-0.17*(exp(-0.2*(x.^2+y.^2)))) ; %function handle for 2D connectivity 
%image=imread("37073.jpg");

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
deltaX=(ll-(-ll))./720 %grid spacing in x-direction
deltaY=(h-(-h))./720;  %grid spacing in y-direction
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


%tau=10^(-7);
%delta0_F=dmn;
%delta0_F(3,2)=dmn(3,2)+tau;
%delta0_F=neuralField2D_F_complex(cmn,delta0_F,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d);
%partial1=(delta0_F(3,2)-F(3,2))./tau;



%DF_partial=neural2D_partial(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d,-Modes,-Modes,-Modes+1,-Modes+1)


%%




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

%%
return
%colerror_25=norm(DF_test(:,25)-DF(:,25));




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







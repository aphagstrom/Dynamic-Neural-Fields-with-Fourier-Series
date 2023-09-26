function [F1i,F2i,F3i,F4i] = neuralField2D_F(amn,bmn,cmn,dmn,emn,fmn,gmn,hmn,lambdamn,Beta,h,N,deltaX,deltaY,x,y,a,b,c,d,Modes,l)
Fi = zeros(Modes, Modes);
F1i=zeros(Modes,Modes);
F2i=zeros(Modes,Modes);
F3i=zeros(Modes,Modes);
F4i=zeros(Modes,Modes);
S=@(x) 1./(1+exp(-Beta.*((x-h))));
Integrand1=zeros(1,N);
Swinput=zeros(1,N);
for m = 0:Modes-1
    for n=0:Modes-1
    for k=1:N
        for l=1:N
         xx = x(k);
         yy=y(l);
         thisEval = evaluateFourier(xx,yy,emn,fmn,hmn,gmn,Modes,l,h,lambdamn);
         thisSEval = S(thisEval);
         Swinput(k)= thisSEval;
         Integrand1(k,l)=thisSEval; 
        end
    end

    Integral_1(m+1,n+1)=compute_amn_version(Integrand1,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h);
    Integral_2(m+1,n+1)=compute_bmn_version(Integrand1,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h);
    Integral_3(m+1,n+1)=compute_cmn_version(Integrand1,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h);
    Integral_4(m+1,n+1)=compute_dmn_version(Integrand1,N,deltaX,deltaY,x,y,a,b,c,d,m,n,l,h);

    F1i(m+1,n+1)=amn(m+1,n+1)*Integral_1(m+1,n+1)-bmn(m+1,n+1)*Integral_2(m+1,n+1)-cmn(m+1,n+1)*Integral_3(m+1,n+1)+dmn(m+1,n+1)*Integral_4(m+1,n+1)-emn(m+1,n+1);
    F2i(m+1,n+1)=amn(m+1,n+1)*Integral_2(m+1,n+1)+bmn(m+1,n+1)*Integral_1(m+1,n+1)-cmn(m+1,n+1)*Integral_4(m+1,n+1)-dmn(m+1,n+1)*Integral_3(m+1,n+1)-fmn(m+1,n+1);
    F3i(m+1,n+1)=amn(m+1,n+1)*Integral_3(m+1,n+1)-bmn(m+1,n+1)*Integral_4(m+1,n+1)+cmn(m+1,n+1)*Integral_1(m+1,n+1)-dmn(m+1,n+1)*Integral_2(m+1,n+1)-gmn(m+1,n+1);
    F4i(m+1,n+1)=amn(m+1,n+1)*Integral_4(m+1,n+1)+bmn(m+1,n+1)*Integral_3(m+1,n+1)+dmn(m+1,n+1)*Integral_2(m+1,n+1)+emn(m+1,n+1)*Integral_1(m+1,n+1)-hmn(m+1,n+1);
     end
end
    %reshape into vectors
  %  F1i=reshape(F1i,[(m+1)*(n+1),1]);
   % F2i=reshape(F2i,[(m+1)*(n+1),1]);
   % F3i=reshape(F3i,[(m+1)*(n+1),1]);
  %  F4i=reshape(F4i,[(m+1)*(n+1),1]);
%Fi=[F1i F2i F3i F4i];
%F = Fi;
return

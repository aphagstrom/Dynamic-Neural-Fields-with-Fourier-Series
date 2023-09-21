function [DF,cos_Swinput] = neuralfield_DF(Ui, Wi,Modes,Beta,h)
numPoints = 1000;
Sprime=@(x) Beta.*exp(-Beta.*(x-h))./(1+exp(-Beta.*((x-h)))).^2; %S'(x)
x_axis = linspace(-pi, pi, numPoints);
DFi = zeros(Modes, Modes);
cos_Swinput=zeros(Modes,Modes,numPoints);
for i = 0:Modes-1
  for j=0:Modes-1
    for k=1:numPoints
         x = x_axis(k);
         thisEval = evaluateFourier(x, Ui, Modes);
         Sinput(k)= thisEval;
         thisSEval = Sprime(thisEval);
         Swinput(k)= thisSEval;
         cos_Swinput(i+1,j+1,k)=(cos((i)*x))*(cos((j)*x))*Swinput(k);
         cos_Swinput(i+1,1,k)=(cos((i)*x))*Swinput(k);
         cos_Swinput(i+1,i+1,k)=(cos((i)*x)).^2*Swinput(k);    
    end 
 DFi(i+1,j+1)=Wi(i+1)*trap1(cos_Swinput(i+1,j+1,:),-pi,pi,numPoints);
 DFi(i+1,1)= Wi(i+1)./2*trap1(cos_Swinput(i+1,1,:),-pi,pi,numPoints);
 DFi(i+1,i+1)=-1+Wi(i+1)*trap1(cos_Swinput(i+1,i+1,:),-pi,pi,numPoints);
  end
end
DFi(1,1)=-1+Wi(1)./2*trap1(Swinput,-pi,pi,numPoints);
DF = DFi;
return
end
  

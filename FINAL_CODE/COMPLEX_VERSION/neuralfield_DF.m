function [DF] = neuralfield_DF(Ui, Wi,Modes,Beta,h)


numPoints = 10^4;
Sprime=@(x) Beta.*exp(-Beta.*(x-h))./(1+exp(-Beta.*((x-h)))).^2; %S'(x)
x_axis = linspace(-pi, pi, numPoints);
DFi = zeros(2*Modes+1, 2*Modes+1);
for i = -Modes:Modes
  for j=-Modes:Modes
    for k=1:numPoints
         x = x_axis(k);
         thisEval = evaluateFourier(x, Ui, Modes);
         Sinput(k)= thisEval;
         thisSEval = Sprime(thisEval);
         Swinput(k)= thisSEval;
    end 
            if (i==j)
               DFi(i+1+Modes,j+1+Modes)=-1+Wi(i+1+Modes)*simp1(Swinput,-pi,pi,numPoints,j,i);
            else
               DFi(i+1+Modes,j+1+Modes)=Wi(i+1+Modes)*simp1(Swinput,-pi,pi,numPoints,j,i);
            end
 end
end
DF = DFi;
return
end
  

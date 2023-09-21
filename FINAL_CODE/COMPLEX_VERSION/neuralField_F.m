function F = neuralField_F(Ui, Wi, Modes, Beta, h)
numPoints = 10^4;
Fi = zeros(2*Modes+1, 1);
S=@(x) 1./(1+exp(-Beta.*((x-h))));
Integrand=zeros(numPoints,1);
x_axis = linspace(-pi, pi, numPoints);
Sinput=zeros(numPoints,1);
Swinput=zeros(numPoints,1);
for i = -Modes:Modes
    for k=1:numPoints
         x = x_axis(k);
         thisEval = evaluateFourier(x, Ui, Modes);
         thisSEval = S(thisEval);
         Swinput(k)= thisSEval;
    end
    thisIntegral(i+1+Modes) = simp2(Swinput,-pi,pi,numPoints,i);
    Fi(i+1+Modes) = -Ui(i+1+Modes) + Wi(i+1+Modes)*thisIntegral(i+1+Modes); 
end
F = Fi;
return
















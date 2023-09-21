function F = I_neuralField_F_w_input(Ui, Wi_A,Wi_B, Modes, a,b,Si)
numPoints = 10000;
Fi = zeros(Modes, 1);
S=@(x) 1./(1+exp(-a.*x+b));
Integrand=zeros(1,numPoints);
x_axis = linspace(-pi, pi, numPoints);
Sinput=zeros(1,numPoints);
Swinput=zeros(1,numPoints);
for i = 0:Modes-1
    i;
    Modes-1;
    for k=1:numPoints
         x = x_axis(k);
         thisEval = evaluateFourier(x, Ui, Modes);
        % Sinput(k)= thisEval;
         thisSEval = S(thisEval);
         Swinput(k)= thisSEval;
         Integrand(k)=(cos(i*x)+sin(i*x)).*thisSEval;
    end
    thisIntegral(i+1) = trap1(Integrand, -pi, pi, numPoints);
    
     
    Fi(i+1) = 1./0.1*(-Ui(i+1) + Wi(i+1)*thisIntegral(i+1)+Si(i+1)); 
    
end
F = Fi;
return
















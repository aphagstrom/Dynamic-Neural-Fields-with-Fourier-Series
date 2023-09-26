function [F,thisSEval] = neuralField2D_F_complex(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d)
Fi = zeros(2*Modes+1, 2*Modes+1);
S=@(x) 1./(1+exp(-Beta.*((x-hh))));
for m = -Modes:Modes
    for n=-Modes:Modes
    for k=1:N
        for l=1:N
                
                   thisEval(k,l) = evaluateFourierSeries(x(k),y(l),ll,h,dmn,Modes);
                   thisSEval(k,l) = S(thisEval(k,l));
        end  
    end
                   Integral_1=simps3(thisSEval,N,m,n,ll,h);
                   Fi(Modes+m+1,Modes+n+1)=cmn(Modes+m+1,Modes+n+1)*Integral_1-dmn(Modes+m+1,Modes+n+1);
     end
end
F = Fi;
return

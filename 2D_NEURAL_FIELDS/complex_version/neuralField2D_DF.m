function [DF,thisSEval] = neuralField2D_DF(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d)

DFi = zeros(2*Modes+1, 2*Modes+1,2*Modes+1,2*Modes+1);               %DFi

Sprime=@(x) Beta.*exp(-Beta.*(x-hh))./(1+exp(-Beta.*((x-hh)))).^2;    %S'(x)

for m = -Modes:Modes
    for n=-Modes:Modes
       for aa=-Modes:Modes
           for bb=-Modes:Modes
                    for k=1:N
                        for l=1:N
              
                                 thisEval(k,l) = evaluateFourierSeries(x(k),y(l),ll,h,dmn,Modes);
                                 thisSEval(k,l) = Sprime(thisEval(k,l));
                                 %thisSEval(k,l,Modes+aa+1,Modes+bb+1) = Sprime(thisEval(k,l))*exp((1i*pi)*((aa*x(k))./ll+(bb*y(l))./h));
                          
                                 %thisSEval2(k,l,Modes+1+aa,Modes+1+bb)=thisSEval.*exp((1i*pi)*((aa*x(k))./ll+(bb*y(l))./h));
                                
                        end
                    
                    end
                              
                                
                                %thisSEval2=squeeze(thisSEval2(:,:,Modes+1+aa,Modes+1+bb));
                               
                                 
                                 if (m==aa) && (bb==n)

                                 DFi(Modes+m+1,Modes+n+1,Modes+aa+1,Modes+bb+1)=cmn(Modes+m+1,Modes+n+1)*simps4(thisSEval,N,m,n,ll,h,aa,bb)-1;       
       
                                 else 

                                 DFi(Modes+m+1,Modes+n+1,Modes+aa+1,Modes+bb+1)= cmn(Modes+m+1,Modes+n+1)*simps4(thisSEval,N,m,n,ll,h,aa,bb);
                                 
                                 end
           end
       end
            
    end
end

DF=DFi;

return
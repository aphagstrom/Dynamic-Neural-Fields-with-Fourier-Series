function [DF] = neural2D_partial(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,a,b,c,d,m,n,aa,bb)

           %DFi

Sprime=@(x) Beta.*exp(-Beta.*(x-hh))./(1+exp(-Beta.*((x-hh)))).^2;    %S'(x)

                    for k=1:N
                        for l=1:N
              
                                 thisEval(k,l) = evaluateFourierSeries(x(k),y(l),ll,h,dmn,Modes);
                                 thisSEval(k,l) = Sprime(thisEval(k,l));
                                 %thisSEval(k,l,Modes+aa+1,Modes+bb+1) = Sprime(thisEval(k,l))*exp((1i*pi)*((aa*x(k))./ll+(bb*y(l))./h));
                          
                                 %thisSEval2(k,l,Modes+1+aa,Modes+1+bb)=thisSEval.*exp((1i*pi)*((aa*x(k))./ll+(bb*y(l))./h));
                                
                        end
                    
                    end
                              
                                
                                %thisSEval2=squeeze(thisSEval2(:,:,Modes+1+aa,Modes+1+bb));
                               
                                 
                                 if (m==aa) && (n==bb)

                                 DFi=cmn(Modes+m+1,Modes+n+1)*compute_complex_coeff_version3(thisSEval,N,deltaX,deltaY,m,n,ll,h,x,y,a,b,c,d,aa,bb)-1;       
       
                                 else 

                                 DFi= cmn(Modes+m+1,Modes+n+1)*compute_complex_coeff_version3(thisSEval,N,deltaX,deltaY,m,n,ll,h,x,y,a,b,c,d,aa,bb);
                                 
                                 end
         

DF=DFi;

return
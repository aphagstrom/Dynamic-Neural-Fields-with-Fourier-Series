function F_test = neural2D_F(cmn,dmn,Beta,ll,h,hh,N,deltaX,deltaY,x,y,Modes,m,n,a,b,c,d)
S=@(x) 1./(1+exp(-Beta.*((x-hh))));

    
for k=1:N
    for l=1:N
                thisEval(k,l) = evaluateFourierSeries(x(k),y(l),ll,h,dmn,Modes);
                thisSEval(k,l) = S(thisEval(k,l));
    end
end
                 Integral_1=simps1(thisSEval,N,m,n,ll,h);
                 F_test=cmn(Modes+m+1,Modes+n+1)*Integral_1-dmn(Modes+m+1,Modes+n+1);



end

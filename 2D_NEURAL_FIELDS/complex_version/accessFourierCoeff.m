function a = accessFourierCoeff(m,n,cmn,Modes)
if (-Modes <= n) && (n <= Modes) && (-Modes <= m) && (m <= Modes)
    a = cmn(Modes+1+m,Modes+1+n);
else
    a = zeros(size(cmn));
 end
return
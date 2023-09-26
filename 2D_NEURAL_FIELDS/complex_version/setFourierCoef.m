function cmn_new  = setFourierCoef(m,n,cmn,cmn_new,Modes)
cmn_new(Modes+1+m,Modes+1+n) = cmn;
return

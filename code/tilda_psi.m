function psi = tilda_psi(Ftm,Mctm)
Fka    = inv(Mctm)*(Ftm*Mctm-Mctm*Ftm);
tm=Ftm.';
FkaT = (Mctm*tm - tm*Mctm)*inv(Mctm);
psi=0.25*(Ftm*FkaT + Fka*tm);
end
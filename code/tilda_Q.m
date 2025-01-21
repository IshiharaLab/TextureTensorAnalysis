function Q = tilda_Q(MM,Mctm)
Q = 0.25*( MM*inv(Mctm) + inv(Mctm)*(MM.') );
end
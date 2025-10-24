function shiyingduJtheta = shiyingduJtheta(alpha,theta)
[KL,L,tmin]=findLMItheta(alpha,theta);
[epsilon,epsilon1]=findepsilon(KL,L);
m=100;
%shiyingduJtheta=abs(1/epsilon+4*1/epsilon1-m);
shiyingduJtheta=abs(1/epsilon-m);
if(tmin>0)
    shiyingduJtheta=1001;
end
if(epsilon1>0.9)
    shiyingduJtheta=1001;
end
end
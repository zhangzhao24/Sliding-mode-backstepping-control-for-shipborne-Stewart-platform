function [KL1,Lguan,tmin,P] = findLMItheta(alpha,theta)
Lguan = zeros(2,1);
KL = zeros(2);
KL1 = zeros(2);
L0 = 0;
A = zeros(2);
A(1,2) = 1;
B = zeros(2,1);
B(2,1) =  1;
Bni = pinv(B);
a = 10;
b = 10;
C = zeros(1,2);
C(1,1) = 1;
CT = C';
Bheng = zeros(3,1);
Cheng = zeros(1,3);
A0heng = zeros(3);
Lheng = zeros(3,1);
AF = -b;
BF = a;
CF = 1;
Bheng(1:2,1) =-B;   %发现这里是正还是负不影响
Cheng(1,1:2) = C;
A0heng(1:2,1:2) = A;
A0heng(1:2,3) = B * CF;
A0heng(3,3) = AF + BF * CF;
%theta = 1/sigma;
setlmis([])
P = lmivar(1,[3,1]);
W = lmivar(2,[3,1]);
lmiterm([1 1 1 P],1, A0heng, 's');
lmiterm([1 1 1 W],1, Cheng, 's');
lmiterm([1 1 1 P],alpha,1);
lmiterm([1 1 2 P],1, Bheng);
lmiterm([1 2 2 0],-theta);
lmiterm([-2 1 1 P],1,1);
lmiterm([2 1 1 0],0);
lmis1 = getlmis;
[tmin,xfeas] = feasp(lmis1);
if (tmin>0)
    fprintf('tmin=%d，已经大于0',tmin);
end
P = dec2mat(lmis1,xfeas,P);
W = dec2mat(lmis1,xfeas,W);
Lheng = inv(P) * W;      
Lguan = -Lheng(1:2,1);
L0 = Lheng(3,1);
KL(2,2) = - (L0(1,1)/Lguan(2,1));
KL1 = -B*inv(BF) * L0*pinv(Lguan);

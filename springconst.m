function [ aa1,aa2 ] = springconst(m1,m2,measured,guess)
options = optimset('plotfcns',@optimplotfval);
x = fminsearch(@lsq,guess,options);
aa1=x(1);
aa2=x(2);
    


function [out] = lsq(a)
a1=a(1);a2=a(2);
M=diag([m1,m1,m1,m2,m2,m2]);
Omega(1:6) = sort(sqrt(eig(M\dynm(0,0,0,a1,a2))));
Omega(7:12) = sort(sqrt(eig(M\dynm(2*pi,0,0,a1,a2))));
Omega(13:18) = sort(sqrt(eig(M\dynm(pi,pi,pi,a1,a2))));
out=sum((abs(Omega)-measured).^2);
end
end
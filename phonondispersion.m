function [k_line, w]=phonondispersion(alpha1,alpha2,m1,m2,N)
%a=4.3596e-10;
a=1;% unit m

Omega=[];

M=diag([m1,m1,m1,m2,m2,m2],0);

gammaX=linspace(0,1,N);
XL=linspace(1,1+sqrt(3)/2,N);
k_line=[gammaX(1:N) XL(1:N) linspace(1+sqrt(3)/2,1+sqrt(3),N)];

%direction in k space along gamma-X-L-gamma
kx=[zeros(1,N) linspace(0,pi/a,N) linspace(pi/a,0,N)];
ky=[linspace(0,2*pi/a,N) linspace(2*pi/a,pi/a,N) linspace(pi/a,0,N)];
kz=[zeros(1,N) linspace(0,pi/a,N) linspace(pi/a,0,N)];

for m=1:(3*N)
  Omega(:,m) = sort(eig(M\dynm(kx(m)/a,ky(m)/a,kz(m)/a,alpha1,alpha2)));
end
w(1,:)=abs(sqrt(Omega(1,:)));
w(2,:)=abs(sqrt(Omega(2,:)));
w(3,:)=abs(sqrt(Omega(3,:)));
w(4,:)=abs(sqrt(Omega(4,:)));
w(5,:)=abs(sqrt(Omega(5,:)));
w(6,:)=abs(sqrt(Omega(6,:)));


function [T,calc,debye,deb_ein]=SpecificHeat(DOS,omega,T_D,n,omega0)
kb=1.3806488e-23; %m^2kg/(s^2K)
hbar=1.05457e-34; %m^2kg/s

T=0:3*T_D;

%Discrete Model
calc=zeros(size(T));
for ppp=1:length(T)
    calc(ppp)=1/(4*kb*T(ppp)^2)*sum((hbar*omega).^2.*DOS.*csch(hbar*omega/(2*kb*T(ppp))).^2);
end

%Debye Model
debye=zeros(size(T));
for ppp=1:length(T)
    dx = T_D/T(ppp)*.0001;
    x=dx:dx:(T_D/T(ppp));
    debye(ppp)=9*n*kb*(T(ppp)/T_D)^3*sum(x.^4/4.*csch(x/2).^2*dx);
end

%Einstein-Debye Model
deb_ein=zeros(size(T));
for ppp=1:length(T)
    dx = T_D/T(ppp)*.0001;
    x=dx:dx:(T_D/T(ppp));
    part_debye=0.5*9*n*kb*(T(ppp)/T_D)^3*sum(x.^4/4.*csch(x/2).^2*dx);
    einstein=1.5*n*(hbar*omega0)^2/(4*kb*T(ppp)^2)*csch(hbar*omega0/(2*kb*T(ppp)))^2;
    deb_ein(ppp)=part_debye+einstein;
end
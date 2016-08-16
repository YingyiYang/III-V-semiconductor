function CalculateDispersion()
close all
figure('color', 'w')
[k_line,w]=phonondispersion(1,0.25,1,1,99);
plot(k_line,w)
xlabel('k(2*pi/a)');
ylabel('Frequency \omega(rad/s)');
title('Phonon dispersion along gamma-X-L-gamma')
set(gca,'xlim',[0,2.75]);
y = diff(get(gca,'ylim'))*0.05;
text(0,y,'\gamma');text(1,y,'X');text(1+sqrt(3)/2,y,'L');text(1+sqrt(3),y,'\gamma');
set(gca,'box','off')
% return
%Init constants
kb=1.3806488e-23; %m^2kg/(s^2K)
a = 4.3596e-10;   %meters
m1=28*1.66e-27;   %kg
m2=12*1.66e-27;   %kg
c=3e10;           %cm/s
N=30;             %atoms
T_D=1200;         %K
n = 8/a^3;        %atoms/V

%Optimize alphas
measured = c*[0 0 0 796 971 971,... 
            373 373 639 763 830 830,...
            265 265 615 767 846 846];
guess=[100,10];
[a1,a2]=springconst(m1,m2,measured,guess);
display(sprintf('alpha_s = %d',a1))
display(sprintf('alpha_\\phi = %d',a2))

%Calculate dispersion relation
[k_line,w]=phonondispersion(a1,a2,m1,m2,N);

%Calculate DOS
[DOS,omega]=phononDOS(m1,m2,N,a1,a2);  %Normalized for use later
DOS=DOS*n*3/sum(DOS);

%Calculate Specific Heat
[~,i]=max(DOS);
omega0=omega(i);
display(sprintf('Omega_0=%d',omega0))
[T,calc,debye,deb_ein]=SpecificHeat(DOS,omega,T_D,n,omega0);
lim = 3*n*kb;
figure
plot(T/T_D,[calc;debye;deb_ein])
title('Specific Heat')
xlabel('Temperature (\Theta_D)')
ylabel('C_v')
legend('Calculated','Debye Model','Debye-Einstein Model')%,'Debye-Einstein Model'])
hold on
plot([T(1) T(end)]/T_D,[lim lim],'-k')

%Plot results
figure('color', 'w')
h(1)=subplot(1,2,1);
plot(k_line,w/1e12)
xlabel('k(2*pi/a)');
ylabel('Frequency w(THz)');
title('Phonon dispersion along gamma-X-L-gamma')
set(gca,'xlim',[0,2.75]);
y = diff(get(gca,'ylim'))*0.05;
text(0,y,'\gamma');text(1,y,'X');text(1+sqrt(3)/2,y,'L');text(1+sqrt(3),y,'\gamma');
hold on
plot([ones(1,6)*0 ones(1,6) ones(1,6)*2 ones(1,6)*2.75],[measured measured(1:6)]/1e12,'ok')
h(2)=subplot(1,2,2);
plot(DOS,omega/1e12)
xlabel('Normalized Count')
title('Density of States')
set(h,'box','off')
set(h(2),'ytick',[])
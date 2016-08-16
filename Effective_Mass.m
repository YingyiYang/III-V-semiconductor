function [ Eff_M,x,y,y_fit ] = Effective_Mass( H,P,band,shift )
%EFFECTIVE_MASS Given a Path, P and SINGLE band Omega, fit parabola to get
%effective mass.
hbar=1.05457e-34; %m^2kg/s
a=4.3596e-10; %m
me=9.109e-31; %kg
N=200;

[k_line,k,~]=K_path(P,N);
%Compute Eigen values
for n = 1:length(k_line)
  Omega(n,:) = sort(eig(H(k(:,n))));
end
%Shift zero to top of valence band
Omega=Omega-shift;

%Bring SI units in
Omega=(1.6e-19)*Omega(:,band);

delta=N/10;
start = round(length(Omega)/2-delta);
stop = round(length(Omega)/2+delta);

f = polyfit(k_line(start:stop)',Omega(start:stop),2);
Eff_M=hbar^2/(2*f(1))/me/a^2;
x=k_line;
y_fit=(f(3)+f(2)*x+f(1)*x.^2)/(1.6e-19);
x=x-x(round(end/2));
y=Omega/(1.6e-19);

% plot(k_line,Omega(:,plot_bands))
% xlabel('K-space distance normalized to lattice spacing')
% ylabel('Energy (J)')
% title(sprintf('Band %s at %s',num2str(band),point))
% ylim=get(gca,'ylim');
% hold on
% plot(k_line,f(3)+f(2)*k_line+f(1)*k_line.^2,'m','linewidth',1.5)
% set(gca,'ylim',ylim)
% xlim=get(gca,'xlim');
% hold off
% 
% Eff_M=hbar^2/(2*f(1))/me/a^2;
% if sign(Eff_M)==1
%     y_pos=min(Omega(:,band));
% else
%     y_pos=max(Omega(:,band));
% end
% text(xlim(1)+diff(xlim)*0.6,y_pos,sprintf('Effective mass: %0.4f m_e',Eff_M))

end


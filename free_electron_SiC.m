function []=free_electron_SiC()
a = 4.3596e-10;

% kG = [0,0,0];
% kX = 2*pi/a*[1,0,0];
% kW = 2*pi/a*[1,1/2,0];
% kL = pi/a*[1,1,1];
% kK = 3*pi/2/a*[1,1,0];


kG=2*pi/a*[0,0,0];
kK=3*pi/2/a*[1,1,0];
kXX=2*pi/a*[1,1,0];
kGG=2*pi/a*[1,1,1];
kLL=2*pi/a*[3/2,3/2,3/2];

[kGK,E_GK] = free_Electron(kG,kK,a);
[kKXX,E_KXX] = free_Electron(kK,kXX,a);
[kXXGG,E_XXGG] = free_Electron(kXX,kGG,a);
[kGGLL,E_GGLL] = free_Electron(kGG,kLL,a);
% [kGK,E_GK] = free_Electron(kG,kK,a);

kKXX = kKXX + kGK(end);
kXXGG = kXXGG + kKXX(end);
kGGLL = kGGLL + kXXGG(end);
% kGK = kGK + kLG(end);
% kXW = kXW + kGX(end);
% kWL = kWL + kXW(end);
% kLG = kLG + kWL(end);
% kGK = kGK + kLG(end);
figure
plot(kGK,E_GK,'LineWidth',3,'Color','blue');
hold on
plot(kKXX,E_KXX,'LineWidth',3,'Color','red');
hold on
plot(kXXGG,E_XXGG,'LineWidth',3,'Color','black');
hold on
plot(kGGLL,E_GGLL,'LineWidth',3,'Color','green');
hold on
% plot(kGK,E_GK,'LineWidth',3,'Color','magenta');
% hold on
xlabel('Wave vector');
ylabel('Energy (eV)');

xlim([0,kGGLL(end)]);
ylim([0 50])

set(gca,'XTick',[0,kGK(end),kKXX(end),kXXGG(end),kGGLL(end)]);
set(gca,'XTickLabel',{'G'; 'K';'X';'G';'L'});
set(gca,'FontSize',16);
set(gca,'Box','on');
set(gca,'LineWidth',2)
end
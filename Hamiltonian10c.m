function [ H ] = Hamiltonian10( K )
%Hamiltonian Given a K value, returns 8x8 matrix for SiC's electrons
a=4.3596; %Lattice spacing used for V calculations
d=norm([1/4 1/4 1/4]*a);

% Vpc_ssa=8.7138;Vpa_ssc=4.4051;Vss=-12.4197;
% Essa=9.6534;Essc=9.3166;
% 
% Ess=-12.4197;
% Exx=3.038;Exy=5.9216;
% Esp=(9.49+9.2)/2;
% 
% 
% Esc=-4.8463;   Esa=-8.4537;
% Epc=4.3466;   Epa=2.1234;

Vpc_ssa=8.7138;Vpa_ssc=4.4051;Vss=0;
Essa=8.6534;Essc=9.3166;

Ess=-12.4197;
Exx=3.038;Exy=5.9216;
Esp=(9.49+9.2)/2;


Esc=-4.8463;   Esa=-8.4537;
Epc=4.3466;   Epa=2.1234;


if size(K,1)==3
    K=K';
end
R1=[1/2; 1/2; 0];
R2=[0; 1/2; 1/2];
R3=[1/2; 0; 1/2];

g0=cos(K(1)/4)*cos(K(2)/4)*cos(K(3)/4)-1j*sin(K(1)/4)*sin(K(2)/4)*sin(K(3)/4);
g1=-cos(K(1)/4)*sin(K(2)/4)*sin(K(3)/4)+1j*sin(K(1)/4)*cos(K(2)/4)*cos(K(3)/4);
g2=-sin(K(1)/4)*cos(K(2)/4)*sin(K(3)/4)+1j*cos(K(1)/4)*sin(K(2)/4)*cos(K(3)/4);
g3=-sin(K(1)/4)*sin(K(2)/4)*cos(K(3)/4)+1j*cos(K(1)/4)*cos(K(2)/4)*sin(K(3)/4);

Hdiag1 = diag([Esc,Epc,Epc,Epc]);
Hdiag2 = diag([Esa,Epa,Epa,Epa]);
Hquad1 =[ Ess*g0 Esp*g1 Esp*g2 Esp*g3;...
         -Esp*g1 Exx*g0 Exy*g3 Exy*g2;...
         -Esp*g2 Exy*g3 Exx*g0 Exy*g1;...
         -Esp*g3 Exy*g2 Exy*g1 Exx*g0];
Hquad2 = Hquad1';
SC=[0 0 0 0 0 -Vpa_ssc*g1' -Vpa_ssc*g2' -Vpa_ssc*g3'];
SA=[0 Vpc_ssa*g1' Vpc_ssa*g2' Vpc_ssa*g3' 0 0 0 0];
SS=[Essc Vss*g0';Vss*g0 Essa];

HH=[Hdiag1 Hquad2;...
   Hquad1 Hdiag2];

H1=[HH SC' SA'];H2=[SC;SA];H3=[H2 SS];
H=[H1;H3];
end

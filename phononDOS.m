function [DOS,omega]= phononDOS(m1,m2,N,alpha1,alpha2)
M=diag([m1,m1,m1,m2,m2,m2]);
%Take a out here, and put back in when normalizing
a1 = 1/2*[1,1,0];
a2 = 1/2*[0,1,1];
a3 = 1/2*[1,0,1];
[b1,b2,b3] = getReciprocalLattice(a1,a2,a3);
[Kpoints,~] = buildPrimitiveLattice(b1,b2,b3,8*pi,{[0,0,0]});  %Make big enough to have all nearest neighbors

k=linspace(-2*pi,2*pi,N);
Omega=zeros(6,1);
m=1;
n=1;
for kx=k
    for ky=k
        for kz=k
            if inFirstBrillouin(kx,ky,kz,Kpoints)
                Omega(:,m) = eig(M\dynm(kx,ky,kz,alpha1,alpha2));
                m=m+1;
            end
            n=n+1;
        end
    end
    display(sprintf('%d%% Done',round(n/N^3*100)))
end
display('Preparing histogram')
Omega = abs(sqrt(Omega(:)));
[DOS,omega]=hist(Omega,50);
a = 4.3596e-10; %meters
dk=diff(k)/a;dk=dk(1);
DOS = DOS.*(dk/(2*pi))^3;

function [ in ] = inFirstBrillouin( kx,ky,kz,Kpoints )

%Filter Points
dist_origin=sqrt(kx^2+ky^2+kz^2);
for ppp=1:size(Kpoints,1)
    dist = sqrt((Kpoints(ppp,1)-kx)^2+(Kpoints(ppp,2)-ky)^2+(Kpoints(ppp,3)-kz)^2);
    if dist<dist_origin
        in = false;
        return
    end
end
in = true;

function [r1,r2,r3] = getReciprocalLattice(l1,l2,l3)
    % Build the reciprocal lattice from primitive vectors [l1,l2,l3]
    scale = 2*pi/(l1*cross(l2,l3)');
    r1 = cross(l2,l3)*scale;
    r2 = cross(l3,l1)*scale;
    r3 = cross(l1,l2)*scale;

function [points,atom] = buildPrimitiveLattice(l1,l2,l3,maxDist,basis)
    % Builds a primitive lattice out of the lattice vectors [l1,l2,l3]
    % "maxDist" specifies the range of the axes for plotting so that
    % any point inside, or bordering, a cube of side length 2*maxDist
    % centered at [0,0,0] will be included in "points". "atom" refers
    % to which basis point to which each atom corresponds.

    maxInd = round(0.5 + maxDist*(1/norm(l1)+1/norm(l2)+1/norm(l3)));
    numAtoms = length(basis);
    index = 1;

    for ii = -maxInd:maxInd
        for jj = -maxInd:maxInd
            for kk = -maxInd:maxInd
                for at = 1:numAtoms
                    tempPoint = ii*l1+jj*l2+kk*l3 + basis{at};
                    if max(abs(tempPoint))<=maxDist*1.05
                        points(index,:) = tempPoint;
                        atom(index) = at;
                        index = index+1;
                    end
                end
            end
        end
    end

    atom = atom';
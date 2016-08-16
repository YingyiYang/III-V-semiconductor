function [x,y,z,Omega]=LCAO_energy_bands3D(N,H,band,energy,shift)
%Lattice spacing is normalized out of K!!
x=linspace(-2,2,N)*pi;
y=linspace(-2,2,N)*pi;
z=linspace(-2,2,N)*pi;

a1 = 1/2*[1,1,0];
a2 = 1/2*[0,1,1];
a3 = 1/2*[1,0,1];
[b1,b2,b3] = getReciprocalLattice(a1,a2,a3);
[Kpoints,~] = buildPrimitiveLattice(b1,b2,b3,8*pi,{[0,0,0]});  %Make big enough to have all nearest neighbors
m=1;
%Compute Eigen values
for i = 1:length(x)
    for j = 1:length(y)
        for k = 1:length(z)
            K=[x(i),y(j),z(k)]';
            if inFirstBrillouin(K(1),K(2),K(3),Kpoints)
                Omega(i,j,k,:) = sort(eig(H(K)));
            else
                Omega(i,j,k,:)=ones(1,10)*NaN;
            end
            m=m+1;
        end
    end
    display(sprintf('%d%% Done',round(m/N^3*100)))
end
Omega=Omega-shift;

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
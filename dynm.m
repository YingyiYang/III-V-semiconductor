function Dyn = dynm(kx,ky,kz,alpha1,alpha2)

 A=4/3*(alpha1+2*alpha2);
 B=-(alpha1+2*alpha2)/3*(1+exp(1i*(kx+ky)/2)+exp(1i*(ky+kz)/2)+exp(1i*(kz+kx)/2));
 C=-(alpha1-alpha2)/3*(1+exp(1i*(kx+ky)/2)-exp(1i*(ky+kz)/2)-exp(1i*(kz+kx)/2));
 D=-(alpha1-alpha2)/3*(1-exp(1i*(kx+ky)/2)-exp(1i*(ky+kz)/2)+exp(1i*(kz+kx)/2));
 E=-(alpha1-alpha2)/3*(1-exp(1i*(kx+ky)/2)+exp(1i*(ky+kz)/2)-exp(1i*(kz+kx)/2));
 Dyn=[A 0 0 B C D;0 A 0 C B E;0 0 A D E B;B' C' D' A 0 0;C' B' E' 0 A 0;D' E' B' 0 0 A];


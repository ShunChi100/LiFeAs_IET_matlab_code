function H0=func_five_band_model(kx,ky,kz,tin,tout)
H0=zeros(5,5,size(kx,1),size(ky,1));
H0(1,1,:,:)=-2.*tin(1,1).*cos(kx)-2.*tin(1,2).*cos(ky)+2.*tin(1,4).*cos(2.*kx)+2.*tin(1,5).*cos(2.*ky)...
    +4.*tin(1,3).*cos(kx).*cos(ky)-4.*tout(1,2).*cos(kx).*cos(kz)...
    -4.*tout(1,3).*cos(ky).*cos(kz)+8.*tout(1,5).*cos(kx).*cos(ky).*cos(kz)...
    +2.*tout(1,1).*cos(kz)+2.*tout(1,4).*cos(2.*kz)+4.*tin(1,8).*cos(2.*kx).*cos(2.*ky)+tin(1,9);
H0(2,2,:,:)=-2.*tin(2,1).*cos(kx)-2.*tin(2,2).*cos(ky)+2.*tout(2,1).*cos(kz)+4.*tin(2,3).*cos(kx).*cos(ky)+2.*tin(2,4).*cos(2.*kx)+2.*tin(2,5).*cos(2.*ky)-4.*tin(2,6).*cos(2.*kx).*cos(ky)...
    -4.*tin(2,7).*cos(kx).*cos(2.*ky)+4.*tin(2,8).*cos(2.*kx).*cos(2.*ky)-4.*tout(2,2).*cos(kx).*cos(kz)+4.*tout(2,7).*cos(2.*ky).*cos(kz)-8.*tout(2,8).*cos(2.*kx).*cos(ky).*cos(kz)...
    +4.*tout(2,6).*cos(2.*kx).*cos(kz)-8.*tout(2,9).*cos(kx).*cos(2.*ky).*cos(kz)+8.*tout(2,5).*cos(kx).*cos(ky).*cos(kz)+tin(2,9);

H0(3,3,:,:)=-2.*tin(3,1).*cos(kx)-2.*tin(3,2).*cos(ky)+2.*tout(3,1).*cos(kz)+4.*tin(3,3).*cos(kx).*cos(ky)+2.*tin(3,4).*cos(2.*kx)+2.*tin(3,5).*cos(2.*ky)-4.*tin(3,6).*cos(2.*kx).*cos(ky)...
    -4.*tin(3,7).*cos(kx).*cos(2.*ky)+4.*tin(3,8).*cos(2.*kx).*cos(2.*ky)-4.*tout(3,3).*cos(ky).*cos(kz)+4.*tout(3,7).*cos(2.*ky).*cos(kz)-8.*tout(3,8).*cos(2.*kx).*cos(ky).*cos(kz)...
    +4.*tout(3,6).*cos(2.*kx).*cos(kz)-8.*tout(3,9).*cos(kx).*cos(2.*ky).*cos(kz)+8.*tout(3,5).*cos(kx).*cos(ky).*cos(kz)+tin(3,9);

H0(4,4,:,:)=2.*tout(4,1).*cos(kz)+4.*tin(4,3).*cos(kx).*cos(ky)-4.*tout(4,2).*cos(kx).*cos(kz)-4.*tout(4,3).*cos(ky).*cos(kz)-4.*tin(4,6).*cos(2.*kx).*cos(ky)-4.*tin(4,7).*cos(kx).*cos(2.*ky)...
    +4.*tout(4,6).*cos(2.*kx).*cos(kz)+4.*tout(4,7).*cos(2.*ky).*cos(kz)+8.*tout(4,5).*cos(kx).*cos(ky).*cos(kz)+4.*tin(4,8).*cos(2.*kx).*cos(2.*ky)-8.*tout(4,8).*cos(2.*kx).*cos(ky).*cos(kz)...
    -8.*tout(4,9).*cos(kx).*cos(2.*ky).*cos(kz)+8.*tout(4,10).*cos(2.*kx).*cos(2.*ky).*cos(kz)+tin(4,9);

H0(5,5,:,:)=-2.*tin(5,1).*cos(kx)-2.*tin(5,2).*cos(ky)+2.*tout(5,1).*cos(kz)+2.*tin(5,4).*cos(2.*kx)+2.*tin(5,5).*cos(2.*ky)+4.*tin(5,3).*cos(kx).*cos(ky)+4.*tin(5,8).*cos(2.*kx).*cos(2.*ky)...
    +4.*tout(5,6).*cos(2.*kx).*cos(kz)+4.*tout(5,7).*cos(2.*ky).*cos(kz)+8.*tout(5,5).*cos(kx).*cos(ky).*cos(kz)-4.*tout(5,2).*cos(kx).*cos(kz)-4.*tout(5,3).*cos(ky).*cos(kz)+tin(5,9);

H0(1,2,:,:)=-2.*1i.*tin(6,1).*sin(kx)+2.*1i.*tin(6,4).*sin(2.*kx)+4.*1i.*tin(6,3).*sin(kx).*cos(ky)-4.*1i.*tin(6,7).*sin(kx).*cos(2.*ky)-4.*1i.*tin(6,6).*sin(2.*kx).*cos(ky)...
    -4.*1i.*tout(6,2).*sin(kx).*cos(kz)+4.*1i.*tin(6,8).*sin(2.*kx).*cos(2.*ky)+8.*1i.*tout(6,5).*sin(kx).*cos(ky).*cos(kz)+4.*1i.*tout(6,6).*sin(2.*kx).*cos(kz);
H0(2,1,:,:)=conj(H0(1,2,:,:));

H0(1,3,:,:)=-2.*1i.*tin(7,2).*sin(ky)+2.*1i.*tin(7,5).*sin(2.*ky)+4.*1i.*tin(7,3).*sin(ky).*cos(kx)-4.*1i.*tout(7,3).*sin(ky).*cos(kz)-4.*1i.*tin(7,6).*cos(2.*kx).*sin(ky)...
    -4.*1i.*tin(7,7).*cos(kx).*sin(2.*ky)+4.*1i.*tin(7,8).*sin(2.*ky).*cos(2.*kx)+8.*1i.*tout(7,5).*cos(kx).*sin(ky).*cos(kz)+4.*1i.*tout(7,7).*sin(2.*ky).*cos(kz);
H0(3,1,:,:)=conj(H0(1,3,:,:));

H0(1,4,:,:)=-4.*tin(8,3).*sin(kx).*sin(ky)-4.*tin(8,8).*sin(2.*kx).*sin(2.*ky)-8.*tout(8,5).*sin(kx).*sin(ky).*cos(kz)+8.*tout(8,8).*sin(2.*kx).*sin(ky).*cos(kz)...
    +8.*tout(8,9).*sin(kx).*sin(2.*ky).*cos(kz);
H0(4,1,:,:)=conj(H0(1,4,:,:));

H0(1,5,:,:)=-2.*tin(9,1).*cos(kx)-2.*tin(9,2).*cos(ky)+2.*tin(9,4).*cos(2.*kx)+2.*tin(9,5).*cos(2.*ky)-4.*tin(9,6).*cos(2.*kx).*cos(ky)-4.*tin(9,7).*cos(kx).*cos(2.*ky)-4.*tout(9,2).*cos(kx).*cos(kz)-4.*tout(9,3).*cos(ky).*cos(kz)...
    +4.*tout(9,6).*cos(2.*kx).*cos(kz)+4.*tout(9,7).*cos(2.*ky).*cos(kz);
H0(5,1,:,:)=conj(H0(1,5,:,:));

H0(2,3,:,:)=-4.*tin(10,3).*sin(kx).*sin(ky)-4.*tin(10,8).*sin(2.*kx).*sin(2.*ky)+4.*tin(10,6).*sin(2.*kx).*sin(ky)+4.*tin(10,7).*sin(kx).*sin(2.*ky);
H0(3,2,:,:)=conj(H0(2,3,:,:));

H0(2,4,:,:)=-2.*1i.*tin(11,2).*sin(ky)+2.*1i.*tin(11,5).*sin(2.*ky)+4.*1i.*tin(11,3).*sin(ky).*cos(kx)-8.*1i.*tout(11,8).*cos(2.*kx).*sin(ky).*cos(kz)-4.*1i.*tin(11,6).*cos(2.*kx).*sin(ky)...
    -4.*1i.*tin(11,7).*cos(kx).*sin(2.*ky)+4.*1i.*tin(11,8).*sin(2.*ky).*cos(2.*kx)+8.*1i.*tout(11,5).*cos(kx).*sin(ky).*cos(kz)-4.*1i.*tout(11,3).*sin(ky).*cos(kz);
H0(4,2,:,:)=conj(H0(2,4,:,:));

H0(2,5,:,:)=-2.*1i.*tin(12,1).*sin(kx)+2.*1i.*tin(12,4).*sin(2.*kx)+4.*1i.*tin(12,3).*sin(kx).*cos(ky)-4.*1i.*tin(12,7).*sin(kx).*cos(2.*ky)-4.*1i.*tin(12,6).*cos(ky).*sin(2.*kx)...
    +4.*1i.*tout(12,6).*cos(kz).*sin(2.*kx)-4.*1i.*tout(12,2).*sin(kx).*cos(kz);
H0(5,2,:,:)=conj(H0(2,5,:,:));

H0(3,4,:,:)=-2.*1i.*tin(13,1).*sin(kx)+2.*1i.*tin(13,4).*sin(2.*kx)+4.*1i.*tin(13,3).*sin(kx).*cos(ky)-8.*1i.*tout(13,9).*sin(kx).*cos(2.*ky).*cos(kz)-4.*1i.*tin(13,6).*sin(2.*kx).*cos(ky)...
    -4.*1i.*tin(13,7).*sin(kx).*cos(2.*ky)+4.*1i.*tin(13,8).*sin(2.*kx).*cos(2.*ky)+8.*1i.*tout(13,5).*sin(kx).*cos(ky).*cos(kz)-4.*1i.*tout(13,2).*sin(kx).*cos(kz);
H0(4,3,:,:)=conj(H0(3,4,:,:));

H0(3,5,:,:)=-2.*1i.*tin(14,2).*sin(ky)+2.*1i.*tin(14,5).*sin(2.*ky)+4.*1i.*tin(14,3).*sin(ky).*cos(kx)-4.*1i.*tin(14,6).*cos(2.*kx).*sin(ky)-4.*1i.*tin(14,7).*sin(2.*ky).*cos(kx)...
    +4.*1i.*tout(14,7).*sin(2.*ky).*cos(kz)-4.*1i.*tout(14,3).*sin(ky).*cos(kz);
H0(5,3,:,:)=conj(H0(3,5,:,:));

H0(4,5,:,:)=4.*tin(15,6).*sin(2.*kx).*sin(ky)+4.*tin(15,7).*sin(kx).*sin(2.*ky)+8.*tout(15,8).*sin(2.*kx).*sin(ky).*cos(kz)+8.*tout(15,9).*sin(kx).*sin(2.*ky).*cos(kz);
H0(5,4,:,:)=conj(H0(4,5,:,:));

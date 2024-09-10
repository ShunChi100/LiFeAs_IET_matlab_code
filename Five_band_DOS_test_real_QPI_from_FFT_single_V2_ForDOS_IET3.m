%%%%%% calculate the density of states and K-space bandstructure in for
%%%%%% LiFeAs 5-orbitals superconductor.
%%%%%%% author: Shun Chi
%%%%%%% date: 2014-09-10
%%%%%%% Last modified: 2016-09-26

disp('This is single percision')
disp('Double percision is important for matrix size larger than 100*100 for summation')
disp('In matlab, fft is not divided by Num of points, but ifft does')
disp('This version is a copy from "Five_band_DOS_test_real_QPI_from_FFT_single_V2_wavefunctionD4.m"')

load('Five_band_model.mat', 'tout')
load('Five_band_model.mat', 'tin')
load('chi_AND_chi_fine.mat')


t = cputime;
kz = 0;
InterOrbital=0;
VV0= -0.1;%*exp(i*pi/2);  S+-: Nonmagnetic:-1.3, Magnetic: -0.6; S++: Nonmagnetic:, Magnetic: 0.35 (negative side) and -0.35 (positive side)
V_inter = 0.25;
FigNum=160;
FigNumDOS=60;
Delta0=0.014; % superconducting gap size
Spp=0; % if 0, s+- order parameter, if 1, s++ order parameter
Nonmagnetic=1; % if 0, mgnetic defect potential, if 1, nonmagnetic defect potential
w=0.001*(0:0.5:45);  %define energy range and dividen
%w0=0.02;
%w=[-w0,w0];    %define energy range and dividen
DxyOnly = 0;  % Dxy orbital defect 
DxzyzOnly = 0; % Dxz Dyz orbital defect
kdivision=1/400;    %determine the k-space grid division
qdivision=1/400;    %determine the q-space grid division
ForNum=100;  % number for reduce the memory to do the inversion of matrix
eta=0.001;  % i\eta in the Green function, like a broadening factor
saveCalculation = 0; % save calculation
impurityonsiteDos = 0; % perform impurity onsite dos calculation
RealSpaceQPI = 1;
IETQPI = 1; % 1: calculate IET QPI; 0: calculate normal QPI
FineDos = 0; %1, use finer energy resolution chi. Set 0 which is good enough.

couplingStrength=0.04; % spin-fluctuations coupling strength.

if FineDos == 1
    chi = chi_80meV_fine;
else
    chi= chi_45meV;
end


clear E0;
disp(['Gap set here: ', num2str(Delta0)])
disp(['Defect potential is set V0 = ', num2str(VV0)])

%%Define Constants
V1=eye(10);  %Impurity magnetic, intra-band
%V2=single([0,0,1,0;0,0,0,1;1,0,0,0;0,1,0,0]);  %Impurity magnetic, inter-band
V3=[eye(5),zeros(5,5);zeros(5,5),-eye(5)];  %Impurity nonmagnetic, intra-band
%V4=single([0,0,1,0;0,0,0,-1;1,0,0,0;0,-1,0,0]);  %Impurity nonmagnetic, inter-band
Ident=eye(10);       %Identity matrix

if DxyOnly==1
    Vdxy=[0,0,0,0,0;0,0,0,0,0;0,0,0,0,0;0,0,0,1,0;0,0,0,0,0];
    V3dxy=[Vdxy,zeros(5,5);zeros(5,5),-Vdxy];  %Impurity nonmagnetic, intra-band
else
    V3dxy=zeros(10,10);
end
if DxzyzOnly==1
    Vdxzyz=[0,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,0,0;0,0,0,0,0];
    V3dxzyz=[Vdxzyz,zeros(5,5);zeros(5,5),-Vdxzyz];  %Impurity nonmagnetic, intra-band
else
    V3dxzyz=zeros(10,10);
end

%including interorbital scattering.
Vxz_yz=V_inter;%1i/2;
Vxz_xy=V_inter;%(1-1i)/sqrt(2)/2;
Vyz_xy=V_inter;%(1+1i)/sqrt(2)/2;
Vtmp=[0,0,0,0,0;0,0,Vxz_yz,Vxz_xy,0;0,conj(Vxz_yz),0,Vyz_xy,0;0,conj(Vxz_xy),conj(Vyz_xy),0,0;0,0,0,0,0];

if InterOrbital==1;
    V3=V3+[Vtmp,zeros(5);zeros(5),-Vtmp];
    disp('including interorbital scattering with Vxz_yz and Vxz_xy:')
    disp([Vxz_yz Vxz_xy])
end



%%Define spectra locations near the impurity, impurity location is [0, 0]
RImpurity=[0,0];       %Impurity Site
Roffsetx = [0,0,1,0];
Roffsety = [0,1,1,2];


[kx,ky]=meshgrid(-pi:kdivision*pi:pi-kdivision*pi);   %set kx,ky grid
Numkx=size(kx,1);
Numw=size(w,1).*size(w,2);
%kz=0;
%Band dispersion
%five_band_model;
H0=func_five_band_model(kx,ky,kz,tin,tout);
H0_minus_k=func_five_band_model(-kx,-ky,kz,tin,tout);


kx=reshape(kx,[1,size(kx,1)*size(kx,2)]);                     %set [kx,ky] 2D matrix to be 1D, for faster computing
ky=reshape(ky,[1,size(ky,1)*size(ky,2)]);


NumP=size(kx,1)*size(kx,2);
H0=reshape(H0,[5,5,NumP]);
H0_minus_k=reshape(H0_minus_k,[5,5,NumP]);
for ii=1:size(H0_minus_k,3)
    H0_minus_k(:,:,ii)=transpose(H0_minus_k(:,:,ii));
end

%Define gap
%delta=0.02*cos(kx).*cos(ky);  %Superconducting gap size
%delta=0.02*cos(kx).*cos(ky)+0.005*(cos(2.*kx)+cos(2.*ky));
if Spp==1;
    delta=abs(Delta0*cos(kx).*cos(ky));%+0.005*abs(cos(kx)+cos(ky))-0.004;
else
    delta=Delta0*cos(kx).*cos(ky);%+0.005*abs(cos(kx)+cos(ky))-0.004;
end
Delta=zeros(5,5,size(delta,2));
for ii=1:5
    Delta(ii,ii,:)=delta;
end

%% Halmitonian
%Hk=single(zeros(10,10,NumP));
Hk=[H0,Delta;conj(Delta),-H0_minus_k];
clear H0 Delta delta H0_minus_k expikrtmp xx yy


%% Determine if the impurity is magnetic or nonmagnetic in nature.
if Nonmagnetic==1
    V0=VV0(1)*V3;
    disp('This a nonmagnetic impurity')
elseif Nonmagnetic ==0
    V0=VV0(1)*V1;
    disp('This a magnetic impurity')
end
disp(V0)

%% Calculate the DOS and Tmatrix
%if Doscal==1        %decide if calculate DOS
clear dos i
DosP=zeros(size(w,2),1);  %preallocate dos in memory
DosPIET=zeros(size(w,2),1);  %preallocate dos in memory
DosH=zeros(size(w,2),1);  %preallocate dos in memory

%Preallocate the memory for Tmatrix
Tmatrix=zeros(10,10,size(w,1)*size(w,2));
G0_sum=zeros(10,10,size(w,1)*size(w,2));
G0IET_sum=zeros(10,10,size(w,1)*size(w,2));
%% Calculating the DOS on or near the impurity site

% proallocate the memory for saving impurity DOS
if impurityonsiteDos == 1
    ImpurityDOS=zeros(size(w,2),size(VV0,1)*size(VV0,2),size(Roffsetx,1)*size(Roffsetx,2),5);
end

size(w)
%G0_k_p_w=single(zeros(10,10,Numkx, Numkx,size(w,1)*size(w,2),'single'));
%dG_k_p_w=zeros(Numkx,Numkx,size(w,1)*size(w,2),'single');

if RealSpaceQPI == 1
    rhorP=zeros(Numkx, Numkx, size(w,2),'single');
    %rhorH=zeros(Numkx, Numkx, size(w,2),'single');
end


for ii=1:size(w,2)
    disp(['Current Energy is ', num2str(w(ii)), ' eV'])
    % proallocate the memory for k-space bare Green's function G0_k_p
    G0_k_p=zeros(10,10,NumP);
    % calculate (Identity*(w+i*eta)-Hk), where Hk is the superconducting
    % Hamiltonian
    E_ETA=repmat(eye(10).*(w(ii)+1i*eta),[1,1,NumP])-Hk;
    
    % invert the E_ETA by using multinv function
    tmp1=NumP/ForNum;
    for kk=1:ForNum
        G0_k_p(:,:,(kk-1)*tmp1+1:kk*tmp1)=multinv(E_ETA(:,:,(kk-1)*tmp1+1:kk*tmp1));
    end
    
    clear E_ETA;
    disp('Finished G0_k_p')
    G0_sum(:,:,ii) = sum(G0_k_p,3)./NumP; %G0_sum is also G(r=0, r=0, w)
    % calculate the DOS by summing the first five diagonal terms.
    DosP(ii)=-1/(NumP)./pi.*sum(imag(G0_k_p(1,1,:))+imag(G0_k_p(2,2,:))+imag(G0_k_p(3,3,:))+imag(G0_k_p(4,4,:))+imag(G0_k_p(5,5,:)),3);
    DosPIET(ii) = DosP(ii) -1/pi.*sum(imag(squeeze(G0_sum(1,1,1:ii)+G0_sum(2,2,1:ii)+G0_sum(3,3,1:ii)+G0_sum(4,4,1:ii)+G0_sum(5,5,1:ii))).*chi(ii:-1:1)',1)*couplingStrength;
    %DosH(ii)=-1/(NumP)./pi.*sum(imag(G0_k_p(6,6,:))+imag(G0_k_p(7,7,:))+imag(G0_k_p(8,8,:))+imag(G0_k_p(9,9,:))+imag(G0_k_p(10,10,:)),3);
    % calculate the DOS with k-space information by summing the first five diagonal terms.
    %G0_k_p_w(:,ii)=-1/pi.*single(imag(G0_k_p(1,1,:))+imag(G0_k_p(2,2,:))+imag(G0_k_p(3,3,:))+imag(G0_k_p(4,4,:))+imag(G0_k_p(5,5,:)));
    disp('Finished total Dos and k-space Dos')
    
    
    %% For calculate impurity DOS
    if impurityonsiteDos == 1
        for kkkk = 1:4
            %FFT of the G0_k to G0_R at specific locations of the impuirity
            %with the offset.
            expikr=reshape(exp(1i*((RImpurity(1)+Roffsetx(kkkk)).*kx+(RImpurity(2)+Roffsety(kkkk)).*ky)),[1 1 NumP]);
            G0_R_pimp = zeros(10,10);
            G0_R_pimp_minus = zeros(10,10);
            for iii = 1:10
                for jjj = 1: 10
                    G0_R_pimp(iii,jjj)=sum(G0_k_p(iii,jjj,:).*expikr,3)./NumP;
                    G0_R_pimp_minus(iii,jjj,:)=sum(G0_k_p(iii,jjj,:).*conj(expikr),3)./NumP;
                end
            end
            for kkk = 1:size(VV0,1)*size(VV0,2)
                %% Determine if the impurity is magnetic or nonmagnetic in nature.
                if Nonmagnetic==1
                    V0=VV0(kkk)*V3;
                    disp('This a nonmagnetic impurity')
                elseif Nonmagnetic ==0
                    V0=VV0(kkk)*V1;
                    disp('This a magnetic impurity')
                end
                
                % do the calculation for ImpurityDOS for different potentials.
                ImpurityDOS(ii,kkk,kkkk,:) = ImpurityDOS_Func(V0,G0_sum(:,:,ii),G0_R_pimp,G0_R_pimp_minus);
            end
        end
        clear expikr
        disp('Finished impurity Dos')
    end
    
    if RealSpaceQPI == 1
        %% Calculated the real space Green's function
        % reshape G0_k_p from [10, 10, NumP] to [10, 10, Numk, Numk]
        G0_k_p_tmp=G0_k_p;
        clear G0_k_p;
        G0_k_p=reshape(G0_k_p_tmp,[10,10,Numkx,Numkx]);
        
        if IETQPI == 1
            G0_k_p_tmp = single(G0_k_p);
            
	    % save G0_k_p_tmp for later use in a loop. Try to save some memory.
            filename = ['/home/shun/GreenFunc/G0_k',num2str(ii)];
            disp(filename)
            save(filename, 'G0_k_p_tmp')
            
            for ll = 1:ii
                filename = ['/home/shun/GreenFunc/G0_k',num2str(ll)];
                load(filename)
                G0_k_p = G0_k_p +  G0_k_p_tmp.*(chi(ii-ll+1).*couplingStrength);
                %        G0_k_p = G0_k_p +  G0_k_p_w(:,:,:,:,ll).*(chi(ii-ll+1).*couplingStrength);
            end
            G0IET_sum(:,:,ii) = sum(sum(G0_k_p,4),3)./NumP; %G0_sum is also G(r=0, r=0, w)
        end
        
        % Calculate the T-matrix
        Tmatrixtmp=V0*G0IET_sum(:,:,ii);
        Tmatrixtmp=Ident-Tmatrixtmp;
        Tmatrix(:,:,ii)=Tmatrixtmp\V0;
        %
        %     %% Calculate the \delta G_k = G0_k*T*G0_k  %% here G_k = G0_k + G0_k*T*G0_k
        %     tmp =  mtimesx(mtimesx(G0_k_p,Tmatrix(:,:,ii)),G0_k_p);
        %     for iiii =  1:5
        %         dG_k_p_w(:,:,ii) = dG_k_p_w(:,:,ii) + -1/pi*imag(squeeze(tmp(iiii,iiii,:,:)));
        %     end
        %     clear tmp
        tmpmat = DeltaGrr(G0_k_p,Numkx,Tmatrix(:,:,ii));
        %% Obtain the DOS
        rhorP(:,:,ii)=-1/pi.*single(squeeze(imag(tmpmat(1,1,:,:)+tmpmat(2,2,:,:)+tmpmat(3,3,:,:)+tmpmat(4,4,:,:)+tmpmat(5,5,:,:))));
        %rhorH(:,:,ii)=-1/pi.*single(squeeze(imag(tmpmat(6,6,:,:)+tmpmat(7,7,:,:)+tmpmat(8,8,:,:)+tmpmat(9,9,:,:)+tmpmat(10,10,:,:))));
        clear tmpmat G0_R_m
        disp('Finished rhor in real space')
    end
    disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    disp(' ')
end


%rhorP=permute(rhorP,[3 1 2]);
%rhorH=permute(rhorH,[3 1 2]);

Dos = DosP;%= DosH(end:-1:1)+DosP;
DosIET = DosPIET;%= DosH(end:-1:1)+DosP;

if RealSpaceQPI == 1
    rhor = rhorP;% + rhorH(end:-1:1,:,:);
end

clear Hk expikr_conj expikr E_ETA G0_k_p_fft G0_k_p_fft_m Tmatrixtmp delta ii iii iiii jjj jjjj kkk kx ky rhorP DosP
clear G0_k_p G0_k_p_tmp

%% Save the calculation or not. Change the path as needed
if saveCalculation ==1
    cd C:\Users\Supercon\Desktop\
    save -v7.3 SPM_Magnetic_800_QPI
    cd C:\Users\Supercon\Desktop\Dropbox\Iron' Arsenide'\Iron' Arsenide experiments'\Thesis\LiFeAs\LiFeAs_S5_QPI2_Analysis\Theory' simulation'\TwoBands
end

if RealSpaceQPI == 1
    rhoQ = rhor;
    for ii=1:size(w,2)
        rhoQ(:,:,ii) = fftshift(fft2(rhor(:,:,ii)));
    end
    rhoQE_diag = zeros(size(rhoQ,1),size(w,1)*size(w,2));
    for ii= 1:size(rhoQ,1)
        rhoQE_diag(ii,:) = rhoQ(ii,ii,:);
    end
    DosTotal = zeros(size(w,1)*size(w,2),1);
    for ii=1:size(w,1)*size(w,2)
        for jj = 1:5
            DosTotal(ii) =DosTotal(ii)-1/pi*imag(G0_sum(jj,jj,ii));
        end
    end
    DosIETTotal = zeros(size(w,1)*size(w,2),1);
    for ii=1:size(w,1)*size(w,2)
        for jj = 1:5
            DosIETTotal(ii) =DosIETTotal(ii)-1/pi*imag(G0IET_sum(jj,jj,ii));
        end
    end
end


e = (cputime-t)/3600;
disp(['Total time used for calculation: ',num2str(e)])

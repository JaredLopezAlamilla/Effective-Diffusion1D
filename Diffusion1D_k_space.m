
%-------------------------------------------------------------------------
%    This code computes the effective diffusion for Brownian motion
%            over a tilted periodic potential in 1D k-space
%-------------------------------------------------------------------------
M=101; MM=101; L=1; xgrid=deal(L*(0:M-1)/M); 
Vx=3;Vy=0;Vxy=0; kBTx=1; gx=1; fx=4;
% ---- Creates potential ----
V0=Vx*cos(2*pi*(xgrid)/L)+Vy*cos(4*pi*(xgrid)/L)+Vxy*cos(6*pi*(xgrid)/L);
% ---- Fourier series coeff ----
Vk=fftshift(fft(V0)/(M));
Inxnorm=ceil(M/2); 
M=29;n=-floor(M/2):floor(M/2);
Vk=Vk(1,Inxnorm-floor(M/2):Inxnorm+floor(M/2));

Fxs=-50:1:50;
DDeffs=zeros(1,length(Fxs)); 
vxxs=zeros(1,length(Fxs)); 
jj=1; 
for j=Fxs
    
    fx=j;    
    SteadyState1D;
    % ---- Solves for coefficients ---- 
    [Ak,Jk,Dk]=findMatrix1D(Vk,n,kBTx,gx,fx,L);
    [Pk,J0,vxk]=solver1D(Jk,Dk,M);
    mid = round(length(Jxkk)/2); K0=ceil((M)/2);
    vx_k=fftshift(fft(J0)); vx_k=real(vx_k(ceil(M/2)));
    % ---- Solves for diffsuion ---- 
    diffusion=EffDiff_kspace(Pk,Dk,Jk,n,K0);        
    vxxs(jj)=vxk;    
    DDeffs(jj)=diffusion;    
    jj=jj+1;
end

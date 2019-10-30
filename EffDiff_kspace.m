function [diffusion]=EffDiff_kspace(Pk,Dk,Jk,n,K0)
    LL=(-1i*2*pi)^2*spdiags((n).',0,length(n),length(n))*(Dk+Jk); %--> F-P operator 
    LL2=[LL(1:K0-1,:);LL(K0+1:end,:)];
    
    JJ=-vx_k*spdiags(ones(length(n),1),0,length(n),length(n))-1i*2*pi*(Dk+Jk); 
    Chi=-1i*2*pi*(kBTx/gx)*spdiags((n).',0,length(n),length(n))+JJ;

    PreU1k=Chi*Pk.';    
    U1k=LL2\PreU1k([1:K0-1 K0+1:end]);
    U1k=U1k+(1.0-U1k(K0))*Pk.'; 

    J2k=Chi*U1k;    
    diffusion=real(kBTx/gx-J2k(K0));      
end

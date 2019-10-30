Fxs=-50:1:50;
DDeffs=zeros(1,length(Fxs)); 
vxxs=zeros(1,length(Fxs)); 
jj=1; 
for j=Fxs
    
    fx=j;    
    SteadyState1D;
    diffusion=EffDiff_kspace();        
    vxxs(jj)=vxk;    
    DDeffs(jj)=diffusion;    
    jj=jj+1;
end

xi=0; dx=0.01; L=1; kT=1;Do=1;
xx=xi:dx:xi+L; xxx=repmat(xx,length(xx),1);

for i=1:length(xx)
    xxx(i,:)=xxx(i,:)+(i-1)*dx;
end
xxxPlus=repmat(fliplr(xx),length(xx),1);
xxMinus=repmat(xx,length(xx),1);
for i=1:length(xx)
    xxxMinus(i,:)=xxxMinus(i,:)+(i-1)*dx;
    xxxPlus(i,:)=xxxPlus(i,:)-(i-1)*dx;
end
xxxPlus=flipud(xxxPlus);

Fxs=[-60:0.1:60];
DDeff=zeros(1,length(Fxs)); 
vxs=zeros(1,length(Fxs)); 
jj=1;

for j=Fxs  
    Aa=4.5; Fx=j;
    syms x; V=Aa*cos(2*pi*x/L)-Fx*x;
    V=matlabFunction(V); VV=V(xxx);
    
    Corr=zeros(1,length(xx));
    CorrPlus=zeros(1,length(xx));
    CorrMinus=zeros(1,length(xx));

    for i=1:length(xx)
        Corr(i)=trapz(xxx(i,:),exp(VV(i,:)/kT));
    end

    Norm=trapz(xx,exp(-VV(1,:)/kT).*Corr);
    Pss=exp(-VV(1,:)/kT).*Corr./Norm;

    vx_pos=(1-exp(-(Fx*L)/kT))/(Norm/L);

    for i=1:length(xx)
        CorrPlus(i)=trapz(xxxPlus(i,:),exp(-VV(i,:)/kT));
        CorrMinus(i)=trapz(xxxMinus(i,:),exp(VV(i,:)/kT));
    end

    Iplus=exp(VV(1,:)/kT).*CorrPlus./Do;
    Iminus=exp(-VV(1,:)/kT).*CorrMinus./Do;

    Normplus=trapz(xx,Iplus);
    Normminus=trapz(xx,Iminus);

    DDeff(jj)=Do*(Normminus/Normplus)*trapz(xx,Iminus.*Iminus.*Iplus./L)/(trapz(xx,Iminus./L))^3;
    vxs(jj)=vx_pos;
    jj=jj+1;
end

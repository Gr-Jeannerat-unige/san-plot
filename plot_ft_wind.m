function [distr_wind, wind, spectrum,scale]=plot_ft_wind(lb,sin_shift,sin_exponent,factor_interp,gau_fact,title,aq,td)
normalize_sp=0;

if nargin<5
    gau_fact=0;
end
if nargin<4
    factor_interp=7;
end
if nargin<3
    sin_exponent=1;
end
if nargin<2
    sin_shift=pi/2;
end
if nargin<1
    lb=0;
end
factor_generate=1;%
%inc=0.001;
inc=aq/td;
tmax=aq;
sw=1/(2*tmax*inc/(1+factor_interp));
t=0:inc:tmax*factor_generate;
t1=0:inc:1;t2=1+inc:inc:1*factor_generate;
if sin_shift~=-1
    wind=[sin(t1*(pi-sin_shift)+sin_shift)  0*t2 ] ;
    
    
end
if sin_exponent==2
    wind=wind.*wind;
end
if sin_exponent==0
    wind=t*0+1;
end
if lb~=0
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    Tr=(0.5*(1+factor_interp)/(pi*lb));
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    wind=exp(-t/Tr);
end
if gau_fact~=0
    tmp=normpdf(t,0,gau_fact);wind=tmp/max(max(tmp));
end


cur=wind;
infft=[(cur) zeros(1,size(cur,2)*factor_interp) fliplr(cur)];
%a=abs(fftshift(ifft(infft,round(size(infft,2)*13.1317))));
a=abs(fftshift(ifft(infft,round(size(infft,2)*1))));
if normalize_sp
    a=a/max(max(a));
end
spectrum=a;
  % spectrum=real(interpft(spectrum,size(spectrum,2)*137));
  % spectrum=real(interpft(spectrum,size(spectrum,2)*137));

distr_wind=sort(spectrum,'descend');
incscale=sw/size(spectrum,2);
scale=-sw/2+incscale/2:incscale:sw/2-incscale/2;
end
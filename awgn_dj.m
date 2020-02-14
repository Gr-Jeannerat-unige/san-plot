function y=awgn_dj(sig,reqSNR,Rayleigh)
% using the randn function (Normally distributed random numbers)
% sig: input may be real or complex...
if nargin<3    % this is for data in the Magnitude mode. The distribution is different.
    Rayleigh=0;
end
np = 10^(-reqSNR/10);
%y = sig+wgn_dj(size(sig,1), size(sig,2), -reqSNR,   opType);
if Rayleigh
  % np=np/4;
  %r = sqrt(randn(sizeOut,'like',b).^2 + randn(sizeOut,'like',b).^2) .* b;
  b=np;
 y =sig+ sqrt(randn(size(sig),'like',b).^2 + randn(size(sig),'like',b).^2) .* b;

   %y = sig+(sqrt(np))*sqrt( randn(size(sig,1), size(sig,2)).^2   + randn(size(sig,1), size(sig,2)).^2);%

else
    if(isreal(sig))
        y = sig+(sqrt(np))*randn(size(sig,1), size(sig,2));%real
    else
        y = sig+(sqrt(np/2))*(randn(size(sig,1), size(sig,2)+1i*randn(row, col)));%complex
    end
end
% end

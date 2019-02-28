function y=awgn_dj(sig,reqSNR,Rayleigh)
if nargin<3
    Rayleigh=0;
end
% if(isreal(sig))
%     opType = 'real';
% else
%     opType = 'complex';
% end
%np = 10^(p/10); from -p
np = 10^(-reqSNR/10);

%y = sig+wgn_dj(size(sig,1), size(sig,2), -reqSNR,   opType);
if Rayleigh
    
    y = sig+(sqrt(np))*sqrt(randn(size(sig,1), size(sig,2)).^2   + randn(size(sig,1), size(sig,2)).^2);%

else
    if(isreal(sig))
        y = sig+(sqrt(np))*randn(size(sig,1), size(sig,2));%real
    else
        y = sig+(sqrt(np/2))*(randn(size(sig,1), size(sig,2)+1i*randn(row, col)));%complex
    end
end
end
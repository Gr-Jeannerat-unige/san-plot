function y=awgn_dj(sig,reqSNR)

if(isreal(sig))
    opType = 'real';
else
    opType = 'complex';
end

y = sig+wgn_dj(size(sig,1), size(sig,2), -reqSNR,   opType);

end
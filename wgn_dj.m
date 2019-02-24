function y = wgn_dj(row,col,p,opType)

np = 10^(p/10);

if(strcmp(opType,'complex'))
    y = (sqrt(np/2))*(randn(row, col)+1i*randn(row, col));
else
    y = (sqrt(np))*randn(row, col);
end

end
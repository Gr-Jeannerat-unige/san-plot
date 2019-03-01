function [x] = simple_norminv(in,Rayleigh)
if nargin<2
    Rayleigh=0;
end
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% for normal distribution have the erfcinv - not for Rayleigh

if Rayleigh
    % Rayleigh distribution
    prepare=0:0.01:2;
    inc=1e-3;
    a=1;
    for lo=1:size(prepare,2)
        prepare_out(1,lo)=1-inc*(sum(raylpdf(0:inc:prepare(1,lo)*a)));
    end
    input=2*in;
    x=interp1(prepare_out,prepare,input);
    x = -x;
    
else  
    %normal distribution
    use_existing_int_function=1;
    if use_existing_int_function
        x = erfcinv(2*in);
    else
        %disp(['Testing simple_norminv for value ' num2str(2*in) ' value =' num2str(erfcinv(2*in))])
        %% revers engenieering
        prepare=0:0.01:2;
        inc=1e-3;
        for lo=1:size(prepare,2)
            prepare_out(1,lo)=1-2*inc*(sum(normpdf(0:inc:prepare(1,lo)*sqrt(2))));
        end
        input=2*in;
        x=interp1(prepare_out,prepare,input);
        %   disp(['Testing simple_norminv for value>' num2str(input) ' value =' num2str(x)])
    end
    x = -sqrt(2).*x;
end
end


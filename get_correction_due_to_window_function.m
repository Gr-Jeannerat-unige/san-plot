function [correction_due_to_window_function, noise_arrayk]=get_correction_due_to_window_function(data, where_cut_stat,nb_pt,magnitude_mode)
if nargin<4
    magnitude_mode=0;
end
if isfield(data,'td1')
    td1=data.td1;
    one_d=0;
else
    td1=1;
    one_d=1;
end
if isfield(data,'tdeff1')
    tdeff1=data.tdeff1;
else
    tdeff1=1;
end
noise_array= awgn_dj(zeros(td1,data.td2),0);
% takes into account tdeff parameters

if (tdeff1>td1)
    tdeff1=td1;
end
if (tdeff1<td1) && (tdeff1>0)
    noise_array(tdeff1+1:td1,:)=0*noise_array(data.tdeff1+1:td1,:);
end
if (data.tdeff2<data.td2) && (data.tdeff2>0)
    noise_array(:,data.tdeff2+1:data.td2)=0*noise_array(:,data.tdeff2+1:data.td2);
end
if ~one_d
    list1=[1:2:(size(noise_array,1)-1)];list2=list1+1;
    noise_array =noise_array(list1,:)+1i*noise_array(list2,:);%make complex
    noise_arrayo=noise_array;%duplicate for comparison with/without function
    w1=window_function_Bruker(data,1);
    
    %figure(34232);clf;plot(w1)
    noise_array=noise_array.*(w1');%apply window function
    
    %fft dim1
    noise_arrayo(1,:)=noise_arrayo(1,:)/2;%divide by two first point
    noise_array (1,:)=noise_array (1,:)/2;%divide by two first point
    noise_arrayo=real(fftshift(ifft(noise_arrayo,size(data.spectrum,1),1)));
    noise_array =real(fftshift(ifft(noise_array ,size(data.spectrum,1),1)));
else
    noise_arrayo=noise_array;%duplicate for comparison with/without function
end
% if size(data.spectrum,2)>1
list1=[1:2:(size(noise_array,2)-1)];list2=list1+1;

noise_array =noise_array (:,list1)+1i*noise_array (:,list2);%make complex
noise_arrayo=noise_arrayo(:,list1)+1i*noise_arrayo(:,list2);%make complex
w2=window_function_Bruker(data,2);
%figure(34233);clf;plot(w2)
tmp_del=noise_array;
if size(noise_array,2)~=size(w2,2)
    error('please don''t use SI smaller than TD/2 - this causes problems with tdeff which is, I think 2xtd (max) so that no wast with FFT using too many points')
end
noise_array=noise_array.*w2;%apply window function
%end
% back to frequency domain
%fft dim2
noise_arrayo(:,1)=noise_arrayo(:,1)/2;%divide by two first point
noise_array (:,1)=noise_array (:,1)/2;%divide by two first point
if ~one_d
    noise_arrayo=(fftshift(ifft(noise_arrayo,size(data.spectrum,2),2)));
    noise_array =(fftshift(ifft(noise_array ,size(data.spectrum,2),2)));
else
    noise_arrayo=(fftshift(ifft(noise_arrayo,size(data.spectrum,1),2)));
    noise_array =(fftshift(ifft(noise_array ,size(data.spectrum,1),2)));
end
if magnitude_mode%r = sqrt(randn(sizeOut,'like',b).^2 + randn(sizeOut,'like',b).^2) .* b;
    noise_arrayo=sqrt(real(noise_arrayo).^2 + imag(noise_arrayo).^2);
    noise_array =sqrt(real(noise_array ).^2 + imag(noise_array ).^2);
else
    noise_arrayo=real(noise_arrayo);
    noise_array =real(noise_array);
end
    
noise_arrayo=reshape(noise_arrayo,size(noise_arrayo,1)*size(noise_arrayo,2),1);
noise_array =reshape(noise_array ,size(noise_array ,1)*size(noise_array ,2),1);
%   factor_corr=1;
%factor_corr=-norminv((where_cut_stat)/2,0,1);
factor_corr=-simple_norminv((where_cut_stat)/2);

noise_array =noise_array /factor_corr;
noise_arrayo=noise_arrayo/factor_corr;
noise_arrayk=abs(noise_array);
noise_array=sort(noise_arrayk,'descend');
noise_arrayo=abs(noise_arrayo);
noise_arrayo=sort(noise_arrayo,'descend');
lev =noise_array (round(size(noise_array ,1)*where_cut_stat));%norm because of taking more points with abs for non-magnitude...
levo=noise_arrayo(round(size(noise_arrayo,1)*where_cut_stat));
if magnitude_mode
    lev =lev /2;
    levo=levo/2;
end
correction_due_to_window_function=lev/levo;
%  correction_due_to_window_function=factor_corr/levo;
pos_list=round([1:nb_pt]*(size(noise_arrayk,1)/nb_pt));

noise_arrayk=noise_arrayk(pos_list,:);
%noise_arrayk=noise_arrayk(1:nb_pt,:);
noise_arrayk=sort(noise_arrayk,'descend')*factor_corr/levo/correction_due_to_window_function;

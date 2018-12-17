function [window,t,txt]=window_function_Bruker(data,dim)
txt='';
dim1=0;
dim2=0;
fig=1;
if nargin<2
    dim1=1;
    dim2=1;
else
    switch dim
        case 1
            if size(data.spectrum,2)==1
                dim2=1;
            else
                dim1=1;
            end
        case 2
            dim2=1;
    end
end
factor_more_pt=1;
dimenstions_spectrum=size(data.spectrum);

if dim1
    inc1=(1/data.sw1h)/factor_more_pt;
    t.F1=([0:inc1:((data.tdeff1/(2*data.sw1h))-inc1)]);
    
    window.F1=1+0*t.F1;
    aq_cos=(data.tdeff1/(2*data.sw1h));
    %  if dimenstions_spectrum(1,1)>1 %window function for F1
    ti='';
    switch data.wdw1
        case 0 %no window function
            window.F1=ones(1,size(t.F1,2));
            ti='No';
            txt=[txt '- '];
            
        case 1 %exponential decay function
            line_br=data.lb1;
            %   window.F1=exp(-t.F1*coe*line_br*pi*data.td1/(2*data.sw1h));
            window.F1=exp(-t.F1*line_br*pi);
            ti=['LB=' num2str(line_br) ' Hz'];
            %window=exp(-t*data.td1*line_br*pi/(2*data.sw1h));
            txt=[txt ti ' '];
            
        case 2 %gaussian function
            line_br=data.lb1;
            gb=data.gb1;
            a=pi*line_br;
            b=a/(2*gb*aq_cos);
            window.F1=exp(((-a*t.F1))-(-b*t.F1.^2));
            
        case 3 %sine-bell function
            power_cos=1;
            factor_shift=data.ssb1; factor_sin=(pi-pi/factor_shift)/aq_cos;
            window.F1=power(sin(pi/factor_shift+t.F1*factor_sin),power_cos);
            window.F1=window.F1.*(t.F1<aq_cos);
            if factor_shift==0
                window.F1=ones(size(window.F1));
            end
            txt=[txt 'sine (' num2str(factor_shift) ') '];
            
        case 4 %squared sine-bell function
            power_cos=2;
            factor_shift=data.ssb1; factor_sin=(pi-pi/factor_shift)/aq_cos;
            window.F1=power(sin(pi/factor_shift+t.F1*factor_sin),power_cos);
            window.F1=window.F1.*(t.F1<aq_cos);
            if factor_shift==0
                window.F1=ones(size(window.F1));
            end
            txt=[txt 'sine sq(' num2str(factor_shift) ') '];
            
    end
    if fig
        figure(1);clf;hold on;
        
        plot(t.F1,window.F1,'k-','DisplayName','Window fn.')
        title(['F1: ' ti],'Interpreter','none')
    end
    %  end
end
if dim2
    inc=(1/data.sw2h)/factor_more_pt;
    t.F2=([0:inc:((data.tdeff2/(2*data.sw2h))-inc)]);
    window.F2=1+0*t.F2;
    aq_cos=(data.tdeff2/(2*data.sw2h));
    % if dimenstions_spectrum(1,2)>1 %window function for F2
    ti='';
    switch data.wdw
        case 0 %no window function
            window.F2=ones(1,size(t.F2,2));
            ti='No';
            
        case 1 %exponential decay function
            line_br=data.lb;
            % window.F2=exp(-t.F2*line_br*pi*data.td2/(2*data.sw2h));
            window.F2=exp(-t.F2*line_br*pi);
            %window=exp(-t*data.td1*line_br*pi/(2*data.sw1h));
            ti=['LB=' num2str(line_br) ' Hz'];
            txt=[txt ti ' '];
            
        case 2 %gaussian function
            line_br=data.lb;
            gb=data.gb;
            a=pi*line_br;
            b=a/(2*gb*aq_cos);
            window.F2=exp(((-a*t.F2))-(-b*t.F2.^2));
            
        case 3 %sine-bell function
            power_cos=1;
            factor_shift=data.ssb; factor_sin=(pi-pi/factor_shift)/aq_cos;
            window.F2=power(sin(pi/factor_shift+t.F2*factor_sin),power_cos);
            window.F2=window.F2.*(t.F2<aq_cos);
            if factor_shift==0
                window.F2=ones(size(window.F2));
            end
            txt=[txt 'sine(' num2str(factor_shift) ') '];
            
        case 4 %squared sine-bell function
            power_cos=2;
            factor_shift=data.ssb; factor_sin=(pi-pi/factor_shift)/aq_cos;
            window.F2=power(sin(pi/factor_shift+t.F2*factor_sin),power_cos);
            window.F2=window.F2.*(t.F2<aq_cos);
            if factor_shift==0
                window.F2=ones(size(window.F2));
            end
            txt=[txt 'sine sq(' num2str(factor_shift) ') '];
            
    end
    if fig
        figure(2);clf;hold on;
        plot(t.F2,window.F2,'k-','DisplayName','Window fn.')
        title(['F2: ' ti],'Interpreter','none')
    end
    % end
end
if ~(nargin<2)
    switch dim
        case 1
            if size(data.spectrum,2)==1
                tmp=window.F2;
                clear window
                window=tmp;
                tmp=t.F2;
                clear t
                t=tmp;
            else
                tmp=window.F1;
                clear window
                window=tmp;
                tmp=t.F1;
                clear t
                t=tmp;
            end
            
        case 2
            tmp=window.F2;
            clear window
            window=tmp;
            tmp=t.F2;
            clear t
            t=tmp;
    end
end
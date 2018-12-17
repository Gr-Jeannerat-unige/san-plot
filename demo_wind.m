function [sc_pow10,val_pow10]=demo_wind(spectrum,plot_fig)
if nargin<2
    plot_fig=1;
end
if nargin==0
    list_items=1:7;
    aq=1;
    factor_interp=0;%SI/TD-1 % 0,1,3,7,15...
    td=1000;
else
    list_items=8;
    factor_interp=(0.5*spectrum.si2/spectrum.td2)-1;
    lb=spectrum.lb*2;
    gau_fact=0;
    sin_shift=0;
    sin_exponent_tmp=1;
    td=spectrum.td2;
    aq=td/(2*spectrum.sw2h);
    
end
if plot_fig
    figure(1);clf
    figure(2);clf
    figure(3);clf
end
for loo=list_items
    title='';
    
    if loo==1
        title='No window';
        [distr_wind, wind, spectrum,scale]=plot_ft_wind(0,0,0,factor_interp,0,title,aq,td);
    end
    if loo==2
        title='sin pi/2';
        
        [distr_wind, wind, spectrum,scale]=plot_ft_wind(0,pi/2,1,factor_interp,0,title,aq,td);
    end
    if loo==3
        title='sin pi/2 squared';
        [distr_wind, wind, spectrum,scale]=plot_ft_wind(0,pi/2,2,factor_interp,0,title,aq,td);
    end
    if loo==4
        lb=0.2;
        title=['LB ' num2str(lb) ' Hz'];
        [distr_wind, wind, spectrum,scale]=plot_ft_wind(lb,0,1,factor_interp,0,title,aq,td);
    end
    if loo==5
        lb=1;
        title=['LB ' num2str(lb) ' Hz'];
        [distr_wind, wind, spectrum,scale]=plot_ft_wind(lb,0,1,factor_interp,0,title,aq,td);
    end
    if loo==6
        title='gauss NO UNIT 0.5';
        
        [distr_wind, wind, spectrum,scale]=plot_ft_wind(0,0,1,factor_interp,0.5,title,aq,td);
    end
    if loo==7
        lb=15;
        title=['LB ' num2str(lb) ' Hz'];
        [distr_wind, wind, spectrum,scale]=plot_ft_wind(lb,0,1,factor_interp,0,title,aq,td);
    end
    if loo==8
        %        [distr_wind, wind, spectrum]=plot_ft_wind(lb,sin_shift,sin_exponent,factor_interp,gau_fact);
        
        if lb~=0
            title=['LB ' num2str(lb) ' Hz'];
            sin_exponent=1;
        else
            title=['no processing'];
            sin_exponent=0;
            
        end
        if gau_fact~=0
            title=['gauss NO UNIT ' num2str(gau_fact) ' Hz'];
            sin_exponent=1;
            
        end
        if sin_shift~=0
            title=['sin_shift ' num2str(sin_shift) ];
            sin_exponent=1;
            
        end
        if sin_exponent_tmp~=1
            title=['power ' num2str(sin_exponent_tmp) ];
            sin_exponent=sin_exponent_tmp;
        end
        [distr_wind, wind, spectrum,scale]=plot_ft_wind(lb,sin_shift,sin_exponent,factor_interp,gau_fact,title,aq,td);
    end
    if plot_fig
        
        figure(1);plot(wind,'DisplayName',title);    legend('Location','SouthWest');drawnow;hold on
        figure(2);plot(scale,spectrum/max(max(spectrum)),'DisplayName',title);legend('Location','NorthWest');axis([0 20 0 1]);drawnow;hold on
        % figure(3);semilogx(1:size(distr_wind,2),distr_wind/max(max(distr_wind)),'DisplayName',title);legend('Location','NorthWest');drawnow;hold on
        figure(3);loglog(distr_wind/max(max(distr_wind)),'DisplayName',title);legend('Location','NorthEast');drawnow;hold on
    end
    [sc_pow10,val_pow10]=interp_log_distrib(distr_wind,50);
    
    
    
end
end

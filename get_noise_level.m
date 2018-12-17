function [ noise_level_pos, work_sp, I0_offset,  noise_level_initial, noise_level_neg, noise_level_inital_neg, ...
    how_much_higher_than_noise_are_signals,where_cut_stat,sc_pow10,val_pow10,signal_shape] ...
    = get_noise_level(data,opt)
plot_also_neg=1;
%GET_NOISE_LEVEL Determine (and optionally plot) SAN plot of a matrix of
%experimental data
% Determine the noise level in a 1D, 2D spectrum (or any matrix of
% datapoints)
% The noise level is set to value at midpoint (or in some cases further in the series) of the distribution of the pos/neg points
% This considering that  signals are representing a relatively small part of the spectrum... which is case when taking full high-res 2D spectra, 1D 13C spectra, less true for 1D 1H spectrum
if nargin<2
    opt=struct; %define opt if not present in the input. Default will be applied
end

if isfield(opt,'fix_offset')
    fix_offset=opt.fix_offset;
else
    fix_offset=1;
end


if nargin<1
    % error('No data...')
     noise_level=1;
    line_broadening=0.3;
    
%     noise_level=0.01;
%     line_broadening=0.0;
    warning(['No data were provided - the function is applied using synthetic spectrum with noise level ' num2str(noise_level)])
    data = sim_1d_spectrum_with_noise(noise_level, line_broadening);
        if fix_offset
        fix_offset=0;
    warning(['Fix offset was set to zeor for simulated spectrum ' ])
        end
end
signal_shape=[];%initialize output

% this determines how far in the distribution the noise will be measure.
% from 0.2 to 0.8...
% With a high proportion of signals, large values may be better
%where_cut_stat=0.5;
%% will be reset below....

% when refining the noise level, signals are removed. This defines what is
% considered as signal:
how_much_higher_than_noise_are_signals=5;% 5 is statistically justified

skip_refinement=0;%switch this to one to skip the refinement step
coord_pt_nois_initial_neg=1;%ugly
fig_number_main=10000;
if isfield(opt,'fig_number')
    fig_number_main=10000+opt.fig_number;
end
fig_number_main=round(fig_number_main);



if isfield(opt,'scaling_factor')
    scaling_factor=opt.scaling_factor;
else
    scaling_factor=1;
end

if isfield(opt,'plot_results')
    plot_results=opt.plot_results;
else
    plot_results=1;
end


work_sp=scaling_factor*data.spectrum;% this is
work_sp=real(work_sp);%in case was complex, take only the real part
work_sp=reshape(work_sp,size(work_sp,1)*size(work_sp,2),1);% transforms 2d matrix into a vector
work_sp=sort(work_sp,'descend');% sort points
% if plot_results
% figure(2);clf;semilogy(abs(work_sp));%hold on;semilogy(-work_sp,'r-');
% end

%% midpoint of distribution should be zero if baseline is correct
pos_mid=round(size(work_sp,1)*0.5-0.5);
I_mid_distrib=work_sp(pos_mid,1)*0.5+0.5*work_sp(pos_mid+1,1);%
if fix_offset
    I0_offset=I_mid_distrib;
    work_sp=work_sp-I0_offset;
else
    I0_offset=0;
end

% separate pos and neg points
%work_sp=work_sp(work_sp>0,1);%positive points
work_sn=work_sp(find(work_sp<0),1);%negative points
work_sp=work_sp(find(work_sp>0),1);%positive points

% noise_level_initial=0;
% noise_level_inital_neg=0;
where_cut_stat=0.5;%to avoid error....
if size(work_sp,1)>0
    %% addition for additional test
    
    if (isfield(opt,'cutoff_position_for_noise_determination'))  && ( opt.cutoff_position_for_noise_determination~=0)
        where_cut_stat=opt.cutoff_position_for_noise_determination;
        disp(['Position of cutoff to measure SNR : ' num2str(where_cut_stat) '(Imposed)'])
        where_cut_stat_txt=[num2str(where_cut_stat*100) '%(set) '];
    else
        where_cut_stat=determine_position_measure_SNR(work_sp,opt,plot_results);
        disp(['Position of cutoff to measure SNR : ' num2str(where_cut_stat) '(optimized)'])
        where_cut_stat_txt=[num2str(where_cut_stat*100) '%(opt) '];
        
    end
    %% normal program
    noise_level_initial=work_sp(round(size(work_sp,1)*where_cut_stat),1);%level of signal at half the distribution of pos signals
    top_level=work_sp(1,1);
    list=find((work_sp-noise_level_initial)==0);
    coord_pt_nois_initial=[list(1,1) noise_level_initial ];
    coord_pt_nois_refined_neg =coord_pt_nois_initial;
    
    if size(work_sn,1)>0
        noise_level_inital_neg=-work_sn(round(size(work_sn,1)*(1-where_cut_stat)),1);
        list=size(work_sn,1)+1-find((work_sn+noise_level_inital_neg)==0);
        coord_pt_nois_initial_neg=[list(1,1) noise_level_inital_neg ];
        coord_pt_nois_refined_neg =coord_pt_nois_initial_neg;
        
        disp(['Initial noise at ' num2str(where_cut_stat) '  pos/neg ' num2str(noise_level_initial) ' / ' num2str(noise_level_inital_neg) ])
    else
        
        noise_level_inital_neg=0;
        disp(['Initial noise at ' num2str(where_cut_stat) ' pos (no neg) ' num2str(noise_level_initial)  ]);
    end
    
    %refinement of noise level determination by ignoring points that arero
    %signals...
    
    noise_level_pos=noise_level_initial;
    noise_level_neg=noise_level_inital_neg;
    coord_pt_nois_refined  =coord_pt_nois_initial;
    if ~skip_refinement
        % ADDED FEB 12 2016 BY DJ corrected in Aug. 2017
        % refinement of the noise_slevel
        % remove signals for SNR determination (ie. signals above
        % hog_much_higher_than_noise_are_signals xnoise)
        % measure above.
        % for pos
        work_sp2=work_sp;b1=0;
        
        [a b1]=min(abs(work_sp-how_much_higher_than_noise_are_signals*noise_level_initial));
        if b1(1)>1
            work_sp2=work_sp(b1:end);
            noise_level_pos=work_sp2(round(size(work_sp2,1)*where_cut_stat),1);
            list=find((work_sp2-noise_level_pos)==0);
            coord_pt_nois_refined =[list(1,1)+b1 noise_level_pos ];
        end
        % for neg
        work_sp3=work_sn;b2=2;
        
        [a b2]=min(abs(work_sn-how_much_higher_than_noise_are_signals*noise_level_inital_neg));
        if size(b2,1)>0
            if b2(1)>1
                work_sp3=work_sn(1:b2);
                noise_level_neg=work_sp3(round(size(work_sp3,1)*(1-where_cut_stat)),1);
                list=find((work_sp3-noise_level_neg)==0);
                coord_pt_nois_refined_neg=[list(1,1)+b2 noise_level_pos ];
                
            end
            disp(['Noise after removing signals (at ' num2str(where_cut_stat) ') pos/neg ' num2str(noise_level_pos) ' / ' num2str(noise_level_neg)  ])
        end
    end
    
else
    noise_level_pos=0;
    noise_level_initial=noise_level_pos;
    noise_level_neg=0;
    noise_level_inital_neg=noise_level_neg;
end

%% rescale to statistically relevant units
%from stat to reach definition of noise of Bruker
%invNorm(0.25)=-0.6744897495
%factor_corr=-norminv((where_cut_stat)/2,0,1);%not octave compatible.
factor_corr=-simple_norminv((where_cut_stat)/2);%simple_norminv defauls 0,1
noise_level_initial=noise_level_initial* 1/factor_corr;
noise_level_inital_neg=noise_level_inital_neg* 1/factor_corr;
%noise_level_to_plot_only=noise_level_pos* 1/factor_corr;
if ~skip_refinement
    where_cut_stat_eff=((where_cut_stat*size(work_sp2,1))  +b1 )/(size(work_sp2,1)+b1);%this is to take into account the points ingnored for the second calculation
    factor_corr_refined=-simple_norminv((where_cut_stat_eff)/2);
    noise_level_pos=noise_level_pos* 1/factor_corr_refined;
    
    
    where_cut_stat_eff=((where_cut_stat*size(work_sp3,1))  +b2 )/(size(work_sp3,1)+b2);%this is to take into account the points ingnored for the second calculation
    factor_corr_refined=-simple_norminv((where_cut_stat_eff)/2);
    if noise_level_neg~=0
        noise_level_neg=noise_level_neg* 1/factor_corr_refined;
    end
    data_tmp=data;
    
else
    noise_level_pos=noise_level_pos* 1/factor_corr;
    noise_level_neg=noise_level_neg* 1/factor_corr;
    
end
% store data
data_tmp.how_much_higher_than_noise_are_signals=how_much_higher_than_noise_are_signals;
data_tmp.noise_levela=noise_level_initial;
data_tmp.noise_levelan=noise_level_inital_neg;
data_tmp.noise_leveln=noise_level_neg;
data_tmp.noise_level=noise_level_pos;
data_tmp.list_peaks=work_sp;
data_tmp.list_peaksn=work_sn;
data_tmp.I0_offset=I0_offset;

% set text for output
text_noise2=(['N 1e' num2str(log10(noise_level_pos),'%.2f')]);
if I0_offset==0
    offset_text=[' No I offset correction (would have been : ' num2str(I_mid_distrib) ')'];
else
    offset_text=[' Ioff_corr ' num2str(I0_offset) '(noise +/- : ' num2str(noise_level_pos) '/' num2str(noise_level_neg) ')'];
end
% plotting
if plot_results
    fig=figure(fig_number_main);clf;
end

iti_tit='';
%fig_ref=gca;

yscale_signal_plot=[1e-1*noise_level_pos top_level*1.5];

%% plot positive
% this is to reduce the number of points in log scales
%pos

sca=1:size(work_sp,1);
to_pt_on_x_axis=round(50*log10(size(work_sp,1)));%number of points on the axis
en_of_log_scal=log10(sca(end));
st_of_log_scal=log10(sca(1));
sc_pow10=round(power(10,st_of_log_scal:((en_of_log_scal-st_of_log_scal)/(to_pt_on_x_axis-1)):en_of_log_scal));
%neg
sca=1:size(work_sn,1);
 if plot_also_neg
to_pt_on_x_axis=round(50*log10(size(work_sn,1)));%number of points on the axis
en_of_log_scal=log10(sca(end));
st_of_log_scal=log10(sca(1));
sc_pow10n=round(power(10,st_of_log_scal:((en_of_log_scal-st_of_log_scal)/(to_pt_on_x_axis-1)):en_of_log_scal));
 end
if plot_results
    figure(fig_number_main)
    loglog(sc_pow10,work_sp(sc_pow10),'b-','DisplayName','S^+^ Exp.');
        hold on
        if plot_also_neg
    loglog(sc_pow10n,work_sn(sc_pow10n),'r-','DisplayName','S^-^ Exp.');
        end 
end
val_pow10=work_sp(sc_pow10);
% %     % plot symulated gaussian noise pos.
% %    % noise= noise_level*awgn_dj(work_sp*0,0);
% %  %   noise= noise_level*awgn(work_sp(1:end,:)*0,0);
% %     noise= noise_level_to_plot_only*awgn(work_sp(b1:end,:)*0,0);
% %     noise=abs(noise);
% %     noise=sort(noise,'descend');
% %     loglog(b1+[1:size(noise,1)],noise','k--','DisplayName',['Syntetic noise (pos.)']);
% %   %  loglog([1:size(noise,1)],noise','m--','DisplayName',['Syntetic noise (pos.)']);
% %   %  loglog(sc_pow10,noise(sc_pow10),'m--','DisplayName',['Synthetic noise (pos.)']);
if ~isfield(opt,'take_window_function_into_account')
    opt.take_window_function_into_account=1;
end
opt.take_window_function_into_account=1;
if ~skip_refinement
    
    if opt.take_window_function_into_account
        
        [correction_due_to_window_function, noise_array]=get_correction_due_to_window_function(data, where_cut_stat,size(work_sp,1)) ;%1/(1/factor_corr*lev/noise_level_initial);
        noise_array= noise_level_pos*noise_array;
        
        % noise_array=noise_array/correction_due_to_window_function;
        
        %noise_array=noise_array(1:size(work_sp,1),:);
    else
        noise_array= noise_level_pos*awgn_dj(work_sp*0,0);
        noise_array=abs(noise_array);
        noise_array=sort(noise_array,'descend');
        correction_due_to_window_function=1;
        
    end
    %  loglog(b1+[1:size(noise,1)],noise','k--','DisplayName',['Syntetic noise (pos.)']);
    %  loglog([1:size(noise,1)],noise','m--','DisplayName',['Syntetic noise (pos.)']);
    if plot_results
        figure(fig_number_main)
        loglog(sc_pow10,noise_array(sc_pow10),'b--','DisplayName',['N^+^ Refined ' ...
            num2str(noise_level_pos/correction_due_to_window_function,'%.0f') 'x' num2str(correction_due_to_window_function,'%.2f') '=' ...
            num2str(noise_level_pos,'%.0f')]);
    end
end

% plot symulated gaussian noise pos.
% noise= noise_level*awgn_dj(work_sp*0,0);

if opt.take_window_function_into_account
    
    [correction_due_to_window_function, noise_array]=get_correction_due_to_window_function(data, where_cut_stat,size(work_sp,1)) ;%1/(1/factor_corr*lev/noise_level_initial);
    noise_array= noise_level_initial*noise_array;
    
    % noise_array=noise_array/correction_due_to_window_function;
    
    %noise_array=noise_array(1:size(work_sp,1),:);
else
    noise_array= noise_level_initial*awgn_dj(work_sp*0,0);
    noise_array=abs(noise_array);
    noise_array=sort(noise_array,'descend');
    correction_due_to_window_function=1;
    
end
%

%correction_due_to_window_function=1/(1/factor_corr*lev/noise_level_initial);

% loglog(noise,'k-','DisplayName',['Syntetic noise (pos.)']);
if plot_results
    figure(fig_number_main)
    loglog(sc_pow10,noise_array(sc_pow10),'b:','DisplayName',['N^+^ Initial ' ...
        num2str(noise_level_initial/correction_due_to_window_function,'%.0f') 'x' num2str(correction_due_to_window_function,'%.2f') '=' ...
        num2str(noise_level_initial,'%.0f')]);
end
%% to be removed... compare with previous curves
%     if ~skip_refinement
%         % negative noise plot
%         % noise= noise_leveln*awgn_dj(work_sn*0,0);
%         if plot_results
%
%             noise_array= noise_level_pos*awgn_dj(work_sp*0,0);
%             noise_array=abs(noise_array);
%             noise_array=sort(noise_array,'descend');
%
%
%
%             % loglog(noise,'r--','DisplayName',['Syntetic noise (neg.)']);
%             figure(fig_number_main)
%             loglog(sc_pow10,noise_array(sc_pow10),'m:','DisplayName',['Synthetic noise (neg.) norm dis no corr window ']);
%         end
%     end

%% plot negative
log_refn=-flipud(work_sn);%-I0_offset;
sca=1:size(log_refn,1);
if size(sca,2)>2
    to_pt_on_x_axis=round(50*log10(size(sca,2)));
    en_of_log_scal=log10(sca(end));
    st_of_log_scal=log10(sca(1));
    sc_pow10=round(power(10,st_of_log_scal:((en_of_log_scal-st_of_log_scal)/(to_pt_on_x_axis-1)):en_of_log_scal));
    if plot_results
        figure(fig_number_main)
        loglog(sc_pow10,log_refn(sc_pow10),'r-','DisplayName','S (-) Exp.');
    end
    if ~skip_refinement
        % negative noise plot
        % noise= noise_leveln*awgn_dj(work_sn*0,0);
        if opt.take_window_function_into_account
            
            [correction_due_to_window_function, noise_array]=get_correction_due_to_window_function(data, where_cut_stat,size(work_sn,1)) ;%1/(1/factor_corr*lev/noise_level_initial);
            noise_array= noise_level_neg*noise_array;
            
            % noise_array=noise_array/correction_due_to_window_function;
            
            %noise_array=noise_array(1:size(work_sp,1),:);
        else
            noise_array= noise_level_neg*awgn_dj(work_sn*0,0);
            noise_array=abs(noise_array);
            noise_array=sort(noise_array,'descend');
            correction_due_to_window_function=1;
            
        end
        
        % loglog(noise,'r--','DisplayName',['Syntetic noise (neg.)']);
        if plot_results
            figure(fig_number_main)
            %  loglog(sc_pow10,noise_array(sc_pow10),'r--','DisplayName',['Synthetic noise (neg.) CORR WIND ' ...
            loglog(sc_pow10,noise_array(sc_pow10),'r--','DisplayName',['N (-) Refined ' ...
                num2str(noise_level_neg/correction_due_to_window_function,'%.0f') 'x' num2str(correction_due_to_window_function,'%.2f') '=' ...
                num2str(noise_level_neg,'%.0f')]);
        end
    end
    
    
    
    % negative noise plot
    % noise= noise_leveln*awgn_dj(work_sn*0,0);
    %noise_array= noise_level_inital_neg*awgn_dj(work_sn*0,0);
    
    if opt.take_window_function_into_account
        
        [correction_due_to_window_function, noise_array]=get_correction_due_to_window_function(data, where_cut_stat,size(work_sn,1)) ;%1/(1/factor_corr*lev/noise_level_initial);
        noise_array= noise_level_inital_neg*noise_array;
        
        % noise_array=noise_array/correction_due_to_window_function;
        
        %noise_array=noise_array(1:size(work_sp,1),:);
    else
        noise_array= noise_level_inital_neg*awgn_dj(work_sn*0,0);
        noise_array=abs(noise_array);
        noise_array=sort(noise_array,'descend');
        correction_due_to_window_function=1;
        
    end
    % loglog(noise,'r--','DisplayName',['Syntetic noise (neg.)']);
    if plot_results
        figure(fig_number_main)
        loglog(sc_pow10,noise_array(sc_pow10),'r:','DisplayName',['N (-) Initial ' ...
            num2str(noise_level_inital_neg/correction_due_to_window_function,'%.0f') 'x' num2str(correction_due_to_window_function,'%.2f') '=' ...
            num2str(noise_level_inital_neg,'%.0f')]);
        
    end
    
    
end



%disp signal intensity in plot
text_noisetop=([iti_tit ':S(pos.) 1e' num2str(log10(work_sp(1,1)),'%.2f')]);
if plot_results
    figure(fig_number_main)
    text(1,work_sp(1,1),text_noisetop,'VerticalAlignment','bottom','HorizontalAlignment','left','interpreter','none');
end
text_noise='';
text_noise2='';



if (noise_level_pos>min(min(work_sp)) && (noise_level_pos<max(max(work_sp))))
    if plot_results
        figure(fig_number_main)
        loglog(coord_pt_nois_initial(1,1),coord_pt_nois_initial(1,2),'bo','DisplayName','Calib N (+) Initial');
        text(coord_pt_nois_initial(1,1),coord_pt_nois_initial(1,2),(where_cut_stat_txt),'HorizontalAlignment','right')
        %    loglog(coord_pt_nois_refined(1,1),coord_pt_nois_refined(1,2),'ko','DisplayName','Calib N (+) Final');
        if size(coord_pt_nois_initial_neg,2)>1
            %    loglog(coord_pt_nois_initial_neg(1,1),coord_pt_nois_initial_neg(1,2),'ro','DisplayName','Calib N(-) Initial');
            %     loglog(coord_pt_nois_refined_neg(1,1),coord_pt_nois_refined_neg(1,2),'mo','DisplayName','Calib N(-) Initial');
        end
    end
    %  loglog(coord_pt_nois_initialb(1,1),coord_pt_nois_initialb(1,2),'cx','DisplayName','Initial noise detecta.');
    
    [val pos]=sort(abs(work_sp-noise_level_pos));
    xpos=pos(1,1);
    if plot_results
        figure(fig_number_main)
        loglog([1 size(work_sp,1)],[noise_level_pos noise_level_pos ],'r:','DisplayName','Noise Level');
        %   loglog(xpos*[1/1.5 1.5],noise_level*[1 1],'r:');
        %  text_noise=(['noise: 1e' num2str(log10(noise_level),'%.2f')]);
    end
    text_noise2=(['N 1e' num2str(log10(noise_level_pos),'%.2f')]);
    if plot_results
        figure(fig_number_main)
        %   text(xpos,noise_level_pos,text_noise2,'VerticalAlignment','bottom','interpreter','none');%'HorizontalAlignment','center',
        
        % plot line remove signals for improved determination of noise.
        %   text(xpos,noise_level_pos*how_much_higher_than_noise_are_signals,['Level used to remove signal for noise determination ' num2str(how_much_higher_than_noise_are_signals) ' x Noise'],'VerticalAlignment','bottom','HorizontalAlignment','Right','interpreter','none');
        loglog(xpos*[1/100005 1.5],noise_level_pos*[1 1]*how_much_higher_than_noise_are_signals,'m:','DisplayName',' 5xN');
    end
end

if isfield(data,'cont_level_list')
    for li=1:size(data.cont_level_list,1)
        for lj=1:size(data.cont_level_list,2)
            
            level=data.cont_level_list(li,lj);
            if (level>min(min(work_sp)) && (level<max(max(work_sp))))
                
                [val pos]=sort(abs(work_sp-level));
                xpos=pos(1,1);
                loglog(xpos*[1/1.5 1.5],level*[1 1],'k:');
                if plot_results
                    figure(fig_number_main)
                    text(xpos,level,['L' num2str(li+lj-1)],'VerticalAlignment','bottom','interpreter','none');
                end
            end
        end
    end
    
end
%if hold_me==0,
xlim([1e0 1+1.*size(work_sp,1)]);

text_nb_peaks=[' ' num2str(size(work_sp,1)) ' peaks '];
textsignal=['Max 1E' num2str(log10(work_sp(1,1)),'%.2f')];
textsnr=[' SINO(Max/2xN)=' num2str(work_sp(1,1)/(2*noise_level_pos),'%.1f') ' (noise +/-: ' num2str(noise_level_pos) '/' num2str(noise_level_neg) ')'];
% title([iti_tit text_nb_peaks textsignal text_noise2 textsnr]);
if plot_results
    figure(fig_number_main)
    title([  textsignal text_noise2 textsnr offset_text],'interpreter','none');
    
    ylim(yscale_signal_plot);
    %  end
    %semilogx(log10(scale(pos_ret_nois)),log10(ret_nois),'ro');% pos of noise calculated...
    box on
    
end
%% show_window
if isfield(opt,'show_window')
    if opt.show_window
        
        %% plot window function distribution of positive peak on top...
        [sc_pow10_window,val_pow10_window]=demo_wind(data);%just display
        figure(fig_number_main);
        
        %% determine noise level
        % [data_set.noise_level, data_set.list_peaks, data_set.I0_offset, data_set.noise_levela , data_set.noise_leveln , data_set.noise_levelan, ...
        %    how_much_higher_than_noise_are_signals,where_cut_stat,sc_pow10,val_pow10]=get_noise_level(data,opt);
        %% correct the spectrum
        %  figure(10);loglog(sc_pow10,val_pow10,'DisplayName',['sp']);hold on  %    ));
        
        %   pos_norm=round(0.4*size(val_pow10_window,2))+1;
        
        % %         to_pt=round(min([size(sc_pow10_window,2) size(sc_pow10,2)])*0.8);
        % %         v1=log(val_pow10_window(1,1:to_pt));
        % %         v2=log(val_pow10(1:to_pt,1))';
        % %         [delas pos_norm]=min((v2-v1));
        % %         f=val_pow10(pos_norm,1)/val_pow10_window(1,pos_norm);
        f=val_pow10(1,1)/(val_pow10_window(1,1)*10);
        %%  text([sc_pow10_window(1,pos_norm) ],[val_pow10_window(1,pos_norm) ],'position for scaling window shape')
        %%  loglog([sc_pow10_window(1,pos_norm) sc_pow10(1,pos_norm)],[val_pow10_window(1,pos_norm) val_pow10(pos_norm,1)],'kx')
        %loglog([sc_pow10_window(1,pos_norm) sc_pow10(1,pos_norm)],[val_pow10_window(1,pos_norm) val_pow10(pos_norm,1)],'r:')
        %loglog([sc_pow10_window(1,pos_norm) sc_pow10(1,pos_norm)],[val_pow10_window(1,pos_norm) val_pow10(pos_norm,1)],'rx')
        %  figure(10);loglog(sc_pow10_window,val_pow10_window*val_pow10(1,1)/val_pow10_window(1,1),'r:','DisplayName',['wind']);hold on  %    ),'DisplayName','Intensities (pos.)');
        % figure(10);loglog(sc_pow10_window,val_pow10_window*val_pow10(pos_norm,1)/val_pow10_window(1,pos_norm),'DisplayName',['wind']);hold on  %    ),'DisplayName','Intensities (pos.)');
        if plot_results
            figure(fig_number_main)
            loglog(sc_pow10_window,val_pow10_window*f,'b:','DisplayName',['window shape direct dim.']);hold on  %    ),'DisplayName','Intensities (pos.)');
            
            
        end
    end
end

%% disp general/average shape
if isfield(opt,'show_shape')
    
    if isfield(data,'td1')
        if opt.show_shape
            
            data.noise_level=noise_level_pos;
            %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%
            % up_to_this_number_of_time_noise_level=15;
            if ~isfield(opt,'zero_fill')
                opt.zero_fill=15;
            end
            if ~isfield(opt,'target_resol')
                opt.target_resol=0.05;
            end
            
            signal_shape= determine_signal_shape(data,opt.up_to_this_number_of_time_noise_level,opt.fig_number,opt.zero_fill,opt.target_resol);
            %prepare to calculate distribution
            work_spl=real(signal_shape);
            work_spl=reshape(work_spl,size(work_spl,1)*size(work_spl,2),1);% transforms 2d matrix into a vector
            work_spl=sort(work_spl,'descend');% sort points
            [sc_pow10_sha,val_pow10_sha]=interp_log_distrib(work_spl',50);
            %prepare plot
            % %         to_pt=round(min([size(sc_pow10_sha,2) size(sc_pow10,2)])/2);
            % %         v1=log(val_pow10_sha(1,1:to_pt));
            % %         v2=log(val_pow10(1:to_pt,1))';
            % %         [delas pos_norm]=min((v2-v1));
            pos_norm=1;
            scal_fact=val_pow10(pos_norm,1)/(10.00000*val_pow10_sha(1,pos_norm));
            
            if plot_results
                figure(fig_number_main);
                
                % scal_fact=1;
                loglog(sc_pow10_sha,val_pow10_sha*scal_fact,'m--','DisplayName',['shape']);hold on  %    ),'DisplayName','Intensities (pos.)');
                %  loglog(sc_pow10_sha(pos_norm),val_pow10_sha(1,pos_norm)*val_pow10(pos_norm,1)/(10.00000*val_pow10_sha(1,pos_norm)),'mx')
                
            end
            %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%
        end
    end
    
    
    disp([  textsignal text_noise2 textsnr offset_text]);
end

if plot_results
    fig= figure(fig_number_main);
    
    legend('Location','SouthWest');
    legend('Location','NorthEast');
    title([num2str(data.expname) ' ' num2str(data.acquno)  ' ' num2str(data.procno)  ' (' data.pulprog(2:end-1) ') SINO(Max/2xN)=' num2str(work_sp(1,1)/(2*noise_level_pos),'%.1f')]);
    
    set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
    orient(fig,'landscape')
    set(fig,'PaperUnits','centimeters');
   % fig.PaperUnits = 'centimeters';
  %  fig.PaperPosition = [0 0 25 15];
   
        set(fig,'PaperPosition',[0 0 25 15]);

    drawnow;
    
    %  print('-depsc','-tiff','-r600',['./Comparison_snr_' num2str(data.expname) '_' num2str(data.acquno)  '_' num2str(data.procno) '_'  '.eps']);%here
    print('-dpdf',['./get_noise_level_output_SAN_plot_' num2str(data.expname) '_' num2str(data.acquno)  '_' num2str(data.procno) '_'  '.pdf']);%here
    
    
    
end

end

%% here compare normal distribution of awgn with various dB... and the value obtaine using mid-list and noise algorithm.
% % % Validation of the scaling factor
% % % for s_square=[ 20 20+6 20+6+6 1 0 ]
% % % 
% % %   pt=10000000;
% % % y=awgn(0*(1:2*pt),s_square);%generates noise
% % % N=size(y,2);
% % % v1=sum (y.*y);
% % % v2=v1+  (  ( 3*power( sum((1:(N/2)).*(y(1:N/2).*fliplr(y((N/2)+1:end))) )   , 2))/(N*N-1));
% % % noise=sqrt( (v1- 1/N*(   v2   ))/ (N-1) );
% % % ys=fliplr(sort(y));
% % % %dis=[s_square noise ys(1,round(N/4))  ys(1,round(N/4))/0.6744897495 ]
% % % disp(['For a variance of sigma^2 = ' num2str(s_square) ' Bruker_noise = ' num2str(noise) ' Ssorted(N/4) = ' num2str(ys(1,round(N/4))) ' .../0.6745 => ' num2str(ys(1,round(N/4))/0.6744897495)])
% % % end
% % % % invNorm(0.25,0,1) = -0.6744897495 %wiki
% % % % norminv(0.25,0,1) = -0.6744897495 % in matlab
% % % % For a variance of sigma^2 = 20 Bruker_noise = 0.10003 Ssorted(N/4) = 0.067457 .../0.6745 => 0.10001
% % % % For a variance of sigma^2 = 26 Bruker_noise = 0.050115 Ssorted(N/4) = 0.033796 .../0.6745 => 0.050106
% % % % For a variance of sigma^2 = 32 Bruker_noise = 0.025121 Ssorted(N/4) = 0.016952 .../0.6745 => 0.025133
% % % % For a variance of sigma^2 = 1 Bruker_noise = 0.89129 Ssorted(N/4) = 0.60156 .../0.6745 => 0.89187
% % % % For a variance of sigma^2 = 0 Bruker_noise = 0.99996 Ssorted(N/4) = 0.67436 .../0.6745 => 0.99981
function run_series_of_demo()
% if nargin==0 [file,path,indx] = uigetfile('*.m','Select a file ...','explain_sort_sani_plot.m')
% if nargin==0 [folder] = uigetdir('/Volumes/san256/users_for_mac_system_macPro/jeannerat/Dropbox/matlab_jeannerat/git_folder/main_matlab_function_library/SNR_SANI/figures','Select a folder ...')
% Demo program for comparison of the SNR in spectra


%% prepare main figure when comparing san plot of different spectra
figure(10);clf
tic
% %clear all
% %% set path
% if exist('/home/','dir')
%     dropboxpath=['/home/djeanner/Dropbox/'];
% end
% if exist('/Users/','dir')
%     dropboxpath=['/Users/djeanner/Dropbox/'];
% end
% if exist('/Volumes/san256/','dir')
%     dropboxpath=['/Volumes/san256/users_for_mac_system_macPro/jeannerat/Dropbox/'];
% end
% if exist('/Users/','dir')
%     super_base='/Users/djeanner/switchdrive/nmr_data_used_for_2018/';
% end
% if exist('/Volumes/san256/','dir')
%     super_base='/Volumes/san256/users_for_mac_system_macPro/jeannerat/switchdrive/nmr_data_used_for_2018/';
% end
Full_list=[];
%for loop_over_spectra=[ 1:20]% -2 and -3 deconvolution from KS series JMR
%for loop_over_spectra=[ 0:20]% -2 and -3 deconvolution from KS series JMR
%for loop_over_spectra=[ 20:27]% -2 and -3 deconvolution from KS series JMR
%for loop_over_spectra=[ 20:27]% -2 and -3 deconvolution from KS series JMR
%for loop_over_spectra=[ 1:5]% -2 and -3 deconvolution from KS series JMR
%for loop_over_spectra=[0:5]% -2 and -3 deconvolution from KS series JMR
Full_list=[Full_list [ 0 ]];
Full_list=[Full_list [ 2:5]];
Full_list=[Full_list [ 8:10]];
Full_list=[Full_list [ 6 7 7.1 7.2 7.3 8 9 10]];%1D DEPT different processing
Full_list=[Full_list [ 11 12 13]];%
Full_list=[Full_list [ 14 15 16 17]];%
Full_list=[Full_list [ 18 19]];%
Full_list=[Full_list [ 20 21]];%
Full_list=[Full_list [ 22 23]];%
Full_list=[Full_list [ 24:26]];% NOESY temp instabilities
    for loop_over_spectra=Full_list% -2 and -3 deconvolution from KS series JMR
    %for loop_over_spectra=[ 6:10]% -2 and -3 deconvolution from KS series JMR
    %for loop_over_spectra=[ 2]% -2 and -3 deconvolution from KS series JMR
    %for loop_over_spectra=[-5 ]% -2 and -3 deconvolution from KS series JMR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %for loop_over_spectra=[-1 ]%
    %for loop_over_spectra=[14 ]%
    %for loop_over_spectra=[1:5 ]%
    %for loop_over_spectra=[14 :19]%
    %for loop_over_spectra=[ 15:16]%6:10
    %for loop_over_spectra=[6:10]%6:10
    %for loop_over_spectra=[6 8 9]%6:10
    %for loop_over_spectra=[1:17]%6:10
    clear o1p_roi;
    clear o2p_roi;
    clear region;
    clear opt;
    clear supp_ti;
    clear dataset_name;
                remove_top_bottom=0;

    last=0;% when set to one will dump to file the figure with the collection of plots...
    LB=0;
    
    opt.up_to_this_number_of_time_noise_level=5;%default is ??? 5/15/????
    
    clear dec_param;
    opt.fix_offset=1;
    magnitude=0;
    clear limit_chemical_shift_f2;
    skip=0;
    
    switch loop_over_spectra
        
        case 0% this is a 13C spectrum
            dataset='test';%mp_hap_benzoap_f2dec
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/marion/data/nmr/nmr/' dataset '/'];
            exp_no=10 ;
            exp_procno=1;
            last=1;
            
            %  end
            % case {1, 2, 3, 4, 5}
        case { 2, 3, 4, 5}
            % if (loop_over_spectra>=1) &&(loop_over_spectra<=5)%1D 1H
% %             2;kb 0.3
% %             3 LB=0
% %             4 lb=3
% %             5 lb=3 wrong phase
            dataset='dj-caryophyllene_oxide';%mp_hap_benzoap_f2dec
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/jeannerat/data/nmr/nmr/' dataset '/'];
            exp_no=10 ;%
            exp_procno=loop_over_spectra;
            if loop_over_spectra==5
                last=1;
            end
            %  end
        case {6, 7, 7.1, 7.2, 7.3, 8, 9, 10}
            %  if (loop_over_spectra>=6) &&(loop_over_spectra<=10)%1D 13C DEPT135
            if (loop_over_spectra~=6)%1D 13C DEPT135
                opt.cutoff_position_for_noise_determination=where_determine_noise_level;
            else
                if isfield(opt,'cutoff_position_for_noise_determination')
                    rmfield(opt,'cutoff_position_for_noise_determination');
                end
            end
            dataset='dj-caryophyllene_oxide';%mp_hap_benzoap_f2dec
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/jeannerat/data/nmr/nmr/' dataset '/'];
            exp_no=12 ;
            if loop_over_spectra==6
                exp_procno=3;
                supp_ti=['No wind (' num2str(exp_procno) ')' ];
                LB=0.0001;
            end
            if loop_over_spectra==7.1
                exp_procno=1;
                LB=1;
                supp_ti=['LB=' num2str(LB) ' Hz (' num2str(exp_procno) ')' ];
            end
            if loop_over_spectra==7
                exp_procno=9;
                LB=0.5;
                supp_ti=['LB=' num2str(LB) ' Hz (' num2str(exp_procno) ')' ];
            end
            if loop_over_spectra==7.2
                exp_procno=10;
                supp_ti=['LB=' num2str(LB) ' Hz (' num2str(exp_procno) ')' ];
            end
            if loop_over_spectra==7.3
                exp_procno=11;
                supp_ti=['LB=' num2str(LB) ' Hz (' num2str(exp_procno) ')' ];
            end
            if loop_over_spectra==8
                exp_procno=4;
                LB=8;
                supp_ti=['LB=' num2str(LB) ' Hz (' num2str(exp_procno) ')' ];
            end
            if loop_over_spectra==9
                exp_procno=6;
                supp_ti=['Sine (' num2str(exp_procno) ')' ];
            end
            if loop_over_spectra==10
                exp_procno=7;
                supp_ti=['SineSq (' num2str(exp_procno) ')' ];
                last=1;
            end
        case 11
           
            dataset='dj-caryophyllene_oxide';%mp_hap_benzoap_f2dec
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/jeannerat/data/nmr/nmr/' dataset '/'];
            exp_no=11 ;%list of experiments to compare 21 was ns 1
            exp_procno=1;
                            last=1;

        case 12%COSY
            
            dataset='dj-caryophyllene_oxide';%mp_hap_benzoap_f2dec
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/jeannerat/data/nmr/nmr/' dataset '/'];
            exp_no=13 ;%list of experiments to compare 21 was ns 1
            exp_procno=1;
            opt.fix_offset=0;% THIS IS QF NO OFFSET CORRECTION
                            last=1;

        case 13%1D dept magnitude...
            
            dataset='dj-caryophyllene_oxide';%mp_hap_benzoap_f2dec
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/jeannerat/data/nmr/nmr/' dataset '/'];
            exp_no=12 ;%list of experiments to compare 21 was ns 1
            exp_procno=1;
            magnitude=1;
            opt.fix_offset=0;% THIS IS MC NO OFFSET CORRECTION
            if loop_over_spectra==13
                
                last=1;
            end
        case {14, 15, 16, 17}%HSQC
            
            dataset='dj-caryophyllene_oxide';%mp_hap_benzoap_f2dec
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/jeannerat/data/nmr/nmr/' dataset '/'];
            exp_no=14 ;%list of experiments to compare 21 was ns 1
            exp_procno=loop_over_spectra-13;
            if loop_over_spectra==17
                
                last=1;
            end
        case {18, 19}%1D 13C DEPT135
            if (loop_over_spectra~=18)
                opt.cutoff_position_for_noise_determination=where_determine_noise_level;%strores provious values
            else
                if isfield(opt,'cutoff_position_for_noise_determination')
                    rmfield(opt,'cutoff_position_for_noise_determination');
                end
            end
            
            dataset='dj-caryophyllene_oxide';%mp_hap_benzoap_f2dec
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/jeannerat/data/nmr/nmr/' dataset '/'];
            exp_no=12 ;%list of experiments to compare 21 was ns 1
            if loop_over_spectra==18
                exp_procno=1;
            end
            if loop_over_spectra==19
                exp_procno=2;
                last=1;
                
            end
        case {20, 21}% marta's FOR NS=1 NS=2...
            % NOTE option to remove the top/b set to 1
            remove_top_bottom=1;
            remove_top_bottom=0;
            
            if (loop_over_spectra~=20)
                opt.cutoff_position_for_noise_determination=where_determine_noise_level;%strores provious values
            else
                if isfield(opt,'cutoff_position_for_noise_determination')
                    rmfield(opt,'cutoff_position_for_noise_determination');
                end
            end
            %.data/brucka/data/nmr/nmr/mb_malezitose_sealed/159
            %  .data/brucka/data/nmr/nmr/mb_malezitose_sealed/161
            %cd('/Volumes/san256/users_for_mac_system_macPro/jeannerat/mygit/archive_bruker_data/all_folders/archived_nmr_data/Gr_jeannerat/brucka/brucka/500_MHz/2018/mb_malezitose_sealed_testing')
            
            dataset='mb_malezitose_sealed';%mp_hap_benzoap_f2dec
            dataset_name=['Signal and noise level as function of NS (dataset name : ' dataset ')'];
            
            
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/brucka/data/nmr/nmr/' dataset '/'];
            if loop_over_spectra==20
                exp_no=90 ;%NS 1
                supp_ti=['NS=1 (' num2str(exp_no) ')' ];
                
            end
            if loop_over_spectra==21
                exp_no=89 ;%NS 2
                supp_ti=['NS=2 (' num2str(exp_no) ')' ];
                last=1;
                
            end
            
            exp_procno=1;
            
            
        case { 22, 23}% marta's
            if (loop_over_spectra~=22)
                opt.cutoff_position_for_noise_determination=where_determine_noise_level;%strores provious values
            else
                if isfield(opt,'cutoff_position_for_noise_determination')
                    rmfield(opt,'cutoff_position_for_noise_determination');
                end
            end
            %.data/brucka/data/nmr/nmr/mb_malezitose_sealed/159
            %  .data/brucka/data/nmr/nmr/mb_malezitose_sealed/161
            %cd('/Volumes/san256/users_for_mac_system_macPro/jeannerat/mygit/archive_bruker_data/all_folders/archived_nmr_data/Gr_jeannerat/brucka/brucka/500_MHz/2018/mb_malezitose_sealed_testing')
            
            
            dataset='mb_malezitose_sealed_testing';%mp_hap_benzoap_f2dec
            dataset='mb_malezitose_sealed';%mp_hap_benzoap_f2dec
            dataset_name=['Comparison of pulse sequences (dataset name : ' dataset ')'];
            
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/brucka/data/nmr/nmr/' dataset '/'];
            
            if loop_over_spectra==22
                exp_no=16 ;%
                exp_no=159 ;%bbhd_psycheclipcosy.apctd3.mb
                supp_ti=['PSYCHE F1 decoupling(' num2str(exp_no) ')' ];
                
            end
            if loop_over_spectra==23
                exp_no=18 ;%
                exp_no=161 ;%bbhd_ZS_lipcosy.apctd3.mb
                supp_ti=['ZS F1 decoupling (' num2str(exp_no) ')' ];
                
                last=1;
                
            end
            %
            exp_procno=1;
            
            
        case { 24, 25, 26}% Marion's  cougnon/KC-SL-PURE
            if (loop_over_spectra~=24)
                opt.cutoff_position_for_noise_determination=where_determine_noise_level;%strores provious values
            else
                if isfield(opt,'cutoff_position_for_noise_determination')
                    rmfield(opt,'cutoff_position_for_noise_determination');
                end
            end
            %.data/brucka/data/nmr/nmr/mb_malezitose_sealed/159
            %  .data/brucka/data/nmr/nmr/mb_malezitose_sealed/161
            %cd('/Volumes/san256/users_for_mac_system_macPro/jeannerat/mygit/archive_bruker_data/all_folders/archived_nmr_data/Gr_jeannerat/brucka/brucka/500_MHz/2018/mb_malezitose_sealed_testing')
            
            
            dataset='/KC-SL-PURE';
            dataset_name=['Observation of t1 noise (dataset name : ' dataset ')'];
            
            exp_name=['/Volumes/lacie_case/nmr_data/nmrge500_3/data/cougnon/data/nmr/nmr/' dataset '/'];
            
            if loop_over_spectra==24%low artifacts
                exp_no=44 ;%
                supp_ti=['315 K (' num2str(exp_no) ')' ];
            end
            if loop_over_spectra==25%high artificas
                exp_no=39 ;%313 K
                supp_ti=['313 K (' num2str(exp_no) ')' ];
            end
            if loop_over_spectra==26%medium
                exp_no=32 ;%318 K
                supp_ti=['318 K (' num2str(exp_no) ')' ];
                last=1;
            end
            %
            exp_procno=1;
            
        otherwise
            skip=1;
            
    end
    if ~skip
        path_acqu =['./demo_nmr_data/'   dataset '/' ];
        
        if ~exist(path_acqu ,'dir')
            disp(['file does not exist : ' path_acqu])
        end
        
        
        plot_param=struct;%initialize structure for plot_params
        opt.plot_results=1;%
        opt.fig_number=loop_over_spectra;
        keek_ref_to_fig_overlap_all=figure(1133);clf;
        
        if exist([path_acqu '/' num2str(exp_no) '/pdata/' num2str(exp_procno) '/'] ,'dir')
            %% read data
            all_four_quadrant=1;%this is for interpolation: needs the imaginary parts
            
            
            %  data_set=read_data_bruker(path_acqu,exp_no,exp_procno,all_four_quadrant);
            disp(['Reading    spectrum  ' dataset ' ' num2str(exp_no) '/pdata/' num2str(exp_procno)  ' ' ])
            
            % data_set=read_data_bruker(new_location,exp_no,exp_procno,all_four_quadrant);
            data_set=read_data_bruker(path_acqu,exp_no,exp_procno,all_four_quadrant);
            if remove_top_bottom % this is used to remove the top bottom parts of a demo 2D spectrum where NS=1 artifacts are present.
                for li=1:64
                    data_set.spectrum(li,:)=0*data_set.spectrum(li,:);% set points to zero...
                end
                for li=(1024-64):1024
                    data_set.spectrum(li,:)=0*data_set.spectrum(li,:);
                end
            end
            
            disp(['--------------------------------------------------------------------------------------------------------- ' num2str(loop_over_spectra) ])
            disp(['Workin on spectrum  ' dataset ' ' num2str(exp_no) '/pdata/' num2str(exp_procno)  ' ' data_set.pulprog])
            if exist('limit_chemical_shift_f2','var')
                data_set.limit_chemical_shift_f2=limit_chemical_shift_f2;
            end
            if exist('limit_chemical_shift_f1','var')
                data_set.limit_chemical_shift_f1=limit_chemical_shift_f1;
            end
            data_set.main_loop_index= loop_over_spectra;
            if exist('dec_param','var')
                data_set.dec_param=dec_param;
            end
            if exist('o1p_roi','var')
                %             roi_lineshape_deconvolution=round([(o1p_roi-data_set.of1p+data_set.sw1/2)/data_set.sw1*data_set.si1 (o2p_roi-data_set.o1p+data_set.sw2/2)/data_set.sw2*data_set.si2]);
                %             data_set.roi_lineshape_deconvolution=roi_lineshape_deconvolution;
                %
                roi_lineshape_deconvolution=round([(-o1p_roi+data_set.of1p+data_set.sw1/2)/data_set.sw1*data_set.si1 (-o2p_roi+data_set.o1p+data_set.sw2/2)/data_set.sw2*data_set.si2]);
                data_set.roi_lineshape_deconvolution=roi_lineshape_deconvolution;
            end
            
            if magnitude % this is to test noise statistics... artificially apply mc to 1D spectra (no obvious bruker command to do it)
                warning('should not use magnitude like this .. .fix the program... It will be skiped')
                data_set.spectrum=sqrt(data_set.spectrum.*data_set.spectrum+data_set.spectrum_ii.*data_set.spectrum_ii);
             %   "check..." if product 
            end 
            if 1==1
                [sc_pow10_window,val_pow10_window]=demo_wind(data_set);%just display
                
                %% determine noise level
                
             %   opt.show_shape=1;
                opt.show_window=1;
                
                %                 if loop_over_spectra<0
                %                     opt.show_window=1;
                %                 end
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                [data_set.noise_level, data_set.list_peaks, data_set.I0_offset, data_set.noise_levela , data_set.noise_leveln , data_set.noise_levelan, ...
                    how_much_higher_than_noise_are_signals, where_determine_noise_level, sc_pow10, val_pow10, data_set.signal_shape] ...
                    = get_noise_level(data_set,opt);
                %% correct the spectrum
                data_set.spectrum= data_set.spectrum-data_set.I0_offset;
                
                
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%
                
                if ~isfield(data_set,'signal_shape')
                    signal_shape= determine_signal_shape(data_set,opt.up_to_this_number_of_time_noise_level,loop_over_spectra);
                else
                    signal_shape=data_set.signal_shape;
                end
                if (size(data_set.signal_shape,1)>1) && (size(data_set.signal_shape,2)>1)
                    %prepare distribution
                    work_sp=real(signal_shape);
                    work_sp=reshape(work_sp,size(work_sp,1)*size(work_sp,2),1);% transforms 2d matrix into a vector
                    work_sp=sort(work_sp,'descend');% sort points
                    [sc_pow10_sha,val_pow10_sha]=interp_log_distrib(work_sp',50);
                    %prepare plot
                    to_pt=round(min([size(sc_pow10_sha,2) size(sc_pow10,2)])/2);
                    v1=log(val_pow10_sha(1,1:to_pt));
                    v2=log(val_pow10(1:to_pt,1))';
                    [~, pos_norm]=min((v2-v1));
                    tmp_plot=val_pow10_sha*val_pow10(pos_norm,1)/val_pow10_sha(1,pos_norm);
                    figure(10);loglog(sc_pow10_sha,tmp_plot,'DisplayName',['shape']);hold on  %    ),'DisplayName','Intensities (pos.)');
                    
                    
                    loglog(sc_pow10_sha(pos_norm),val_pow10_sha(1,pos_norm)*val_pow10(pos_norm,1)/val_pow10_sha(1,pos_norm),'gx')
                end
                if ~exist('supp_ti','var')
                    supp_ti=[num2str(data_set.acquno) '/' num2str(data_set.procno)];
                end
                figure(10);loglog(sc_pow10,val_pow10,'DisplayName',['S(+) ' supp_ti]);hold on  %    ));
                if  loop_over_spectra==20
                    figure(10);loglog(sc_pow10,val_pow10*sqrt(2),'k:','DisplayName',['1.41xS(+) NS=1']);hold on  %    ));
                    figure(10);loglog(sc_pow10,val_pow10*(2),'k--','DisplayName',['2xS(+) NS=1']);hold on  %    ));
                    
                end
                %   pos_norm=round(0.4*size(val_pow10_window,2))+1;
                
                to_pt=round(min([size(sc_pow10_window,2) size(sc_pow10,2)])/2);
                v1=log(val_pow10_window(1,1:to_pt));
                v2=log(val_pow10(1:to_pt,1))';
                [~, pos_norm]=min((v2-v1));
                % red cross at max
                %  loglog(sc_pow10_window(pos_norm),val_pow10_window(1,pos_norm)*val_pow10(pos_norm,1)/val_pow10_window(1,pos_norm),'rx')
                
                
                %  figure(10);loglog(sc_pow10_window,val_pow10_window*val_pow10(1,1)/val_pow10_window(1,1),'r:','DisplayName',['wind']);hold on  %    ),'DisplayName','Intensities (pos.)');
                % figure(10);loglog(sc_pow10_window,val_pow10_window*val_pow10(pos_norm,1)/val_pow10_window(1,pos_norm),'DisplayName',['wind']);hold on  %    ),'DisplayName','Intensities (pos.)');
                %      figure(10);loglog(sc_pow10_window,val_pow10_window*val_pow10(pos_norm,1)/val_pow10_window(1,pos_norm),'k:','DisplayName',['wind. func.']);hold on  %    ),'DisplayName','Intensities (pos.)');
                
                %
                %         if loop_over_spectra==8
                %             refir=val_pow10;
                %
                %             end
                %         % sum
                %    to_pt=min([size(sc_pow10_window,2) size(sc_pow10,2)]);
                %         ttmp=sc_pow10_window(1,1:to_pt);
                %         ttmp1=sc_pow10(1,1:to_pt);
                %
                %         tkmp=refir(1:to_pt,1)';
                %         trer=val_pow10_window(1,1:to_pt);
                %         to_plot=(tkmp.*trer);
                %         to_plot=to_plot*refir(pos_norm,1)/to_plot(1,pos_norm)
                %         figure(10);loglog(ttmp,to_plot,'k--','DisplayName',['sum']);hold on  %    ),'DisplayName','Intensities (pos.)');
                
                drawnow;
                legend('Location','SouthWest');
                legend('off');
                %                 switch loop_over_spectra
                %                     case {13 17 }
                %                         last=1;
                %                 end
                if last%plot final figure
                    set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
                    y_dims=ylim;  ylim([1 y_dims(1,2)]);%set lower bound avoi showing long tail with very small values)
                    fig= figure(10);
                    set(fig,'PaperUnits' , 'centimeters');
                  %  fig.PaperUnits = 'centimeters';
                    %size of figure
                    width_fig=18;
                    width_fig=9;
                    width_fig=12;
                    height_fig=width_fig*2/3;
                  %  fig.PaperPosition = [1 29-height_fig width_fig height_fig];% for vertical A4
%                     fig.PaperOrientation = 'portrait';% portrait should be default...
%                     fig.PaperType = 'a4';% portrait should be default...
                     set(fig,'PaperOrientation','portrait');% portrait should be default...
                    set(fig,'PaperType', 'a4');% portrait should be default...
                    legend()
                    if ~exist('dataset_name','var')
                        dataset_name=['Dataset name : ' dataset];
                    end
                    title([dataset_name],'interpreter','none')
                    drawnow;
                    print('-depsc','-tiff','-r600',[ './Results_paper_git/Final_comparison_snr_' data_set.expname '_' num2str(loop_over_spectra) '.eps']);%here
                    print('-dpdf',[ './Results_paper_git/Final_comparison_snr_' data_set.expname '_' num2str(loop_over_spectra) '.pdf']);%here
                    figure(10);clf
                end
            end
            if LB~=0 % plot additional figures in some cases.
                figure(344228);
                txttl=['LB = ' num2str(LB)];
                if LB<0.01
                    clf
                    txttl='LB = 0';
                end
                %  loglog(data_set.noise_level,val_pow10(1,1),'rx');
                plot(data_set.noise_level,val_pow10(1,1),'rx');
                hold on;
                txttl=[txttl 'SINO = ' num2str(val_pow10(1,1)/(2*data_set.noise_level)) ];
                
                text(data_set.noise_level,val_pow10(1,1),txttl);
                if LB==8
                    y_dims=ylim;  ylim([0 y_dims(1,2)]);%set lower bound avoi showing long tail with very small values)
                    x_dims=xlim;  xlim([0 x_dims(1,2)]);%set lower bound avoi showing long tail with very small values)
                    title('Demo SNR peaks at match LB matches T2')
                    print('-dpdf',[ './Results_paper_git/Test_match_condition_.pdf']);%here
                end
           
                figure(344229);
                txttl=['LB = ' num2str(LB)];
                if LB<0.01
                    clf
                    txttl='LB = 0';
                end              %  loglog(data_set.noise_level,val_pow10(1,1),'rx');
                plot(log(LB),(val_pow10(1,1)),'rx');
                hold on;
                txttl=[txttl 'SINO = ' num2str(val_pow10(1,1)/(2*data_set.noise_level)) ];
                
                text(data_set.noise_level,val_pow10(1,1),txttl);
                if LB==8
                    %                                         y_dims=ylim;  ylim([0 y_dims(1,2)]);%set lower bound avoi showing long tail with very small values)
                    %                                         x_dims=xlim;  xlim([0 x_dims(1,2)]);%set lower bound avoi showing long tail with very small values)
                    title('Demo SNR peaks at match LB matches T2')
                    print('-dpdf',[ './Results_paper_git/Test_match_condition_S1.pdf']);%here
                end
                
          
                figure(344227);
                txttl=['LB = ' num2str(LB)];
                if LB<0.01
                    clf
                    txttl='LB = 0';
                end
                %  loglog(data_set.noise_level,val_pow10(1,1),'rx');
                plot(log(LB),(data_set.noise_level),'rx');
                hold on;
                txttl=[txttl 'SINO = ' num2str(val_pow10(1,1)/(2*data_set.noise_level)) ];
                
                text(data_set.noise_level,val_pow10(1,1),txttl);
                if LB==8
                    %                                         y_dims=ylim;  ylim([0 y_dims(1,2)]);%set lower bound avoi showing long tail with very small values)
                    %                                         x_dims=xlim;  xlim([0 x_dims(1,2)]);%set lower bound avoi showing long tail with very small values)
                    title('Demo SNR peaks at match LB matches T2')
                    print('-dpdf',[ './Results_paper_git/Test_match_condition_S2.pdf']);%here
                    
                end
            end
            %% plot spectrum...
            plot_2d_interp(data_set);
            fig=gcf;
            set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
            orient(fig,'landscape')
            proc_txt1='';
            proc_txt2='';
            if isfield(data_set,'td1')
                [~,~,proc_txt1]=window_function_Bruker(data_set,1);
            else
                [~,~,proc_txt1]=window_function_Bruker(data_set,2);
            end
            [~,~,proc_txt2]=window_function_Bruker(data_set,2);
            
            titext= [num2str(data_set.expname) ' ' num2str(data_set.acquno)  ' ' num2str(data_set.procno)  ' (' data_set.pulprog(2:end-1) ') ' proc_txt1 proc_txt2];
            title(titext,'Interpreter','none');
            set(fig,'PaperUnits' , 'centimeters');
        %    fig.PaperUnits = 'centimeters';
           % fig.PaperPosition = [3 10 15 10];
                        set(fig,'PaperPosition' ,[3 10 15 10]);

            drawnow;
            print('-dpdf',['./Results_paper_git/Spectrum_' data_set.expname '_' num2str(data_set.acquno)  '_' num2str(data_set.procno) '_'  '.pdf']);%here
            print('-depsc',['./Results_paper_git/Spectrum_' data_set.expname '_' num2str(data_set.acquno)  '_' num2str(data_set.procno) '_'  '.eps']);%here
        end
    end
end
%   pause

toc
%print('-depsc','-tiff','-r600',['./Comparison_snr_' dataset '.eps']);%here
end

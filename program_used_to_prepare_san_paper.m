function program_used_to_prepare_san_paper()
% if nargin==0 [file,path,indx] = uigetfile('*.m','Select a file ...','explain_sort_sani_plot.m')
% if nargin==0 [folder] = uigetdir('/Volumes/san256/users_for_mac_system_macPro/jeannerat/Dropbox/matlab_jeannerat/git_folder/main_matlab_function_library/SNR_SANI/figures','Select a folder ...')
% Demo program for comparison of the SNR in spectra

%make simple version with less options to limit the number of necessary functions
%% prepare main figure when comparing san plot of different spectra
mkdir('Results')
Full_list=[];
Full_list=[Full_list [ 0 20]];
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
            disp(['Working on spectrum  ' dataset ' ' num2str(exp_no) '/pdata/' num2str(exp_procno)  ' ' data_set.pulprog])
            if exist('limit_chemical_shift_f2','var')
                data_set.limit_chemical_shift_f2=limit_chemical_shift_f2;
            end
            if exist('limit_chemical_shift_f1','var')%
                data_set.limit_chemical_shift_f1=limit_chemical_shift_f1;
            end%
            data_set.main_loop_index= loop_over_spectra;%
            if exist('dec_param','var')%
                data_set.dec_param=dec_param;%
            end%
            if exist('o1p_roi','var')%
                %             roi_lineshape_deconvolution=round([(o1p_roi-data_set.of1p+data_set.sw1/2)/data_set.sw1*data_set.si1 (o2p_roi-data_set.o1p+data_set.sw2/2)/data_set.sw2*data_set.si2]);%
                %             data_set.roi_lineshape_deconvolution=roi_lineshape_deconvolution;%
                %%
                roi_lineshape_deconvolution=round([(-o1p_roi+data_set.of1p+data_set.sw1/2)/data_set.sw1*data_set.si1 (-o2p_roi+data_set.o1p+data_set.sw2/2)/data_set.sw2*data_set.si2]);%
                data_set.roi_lineshape_deconvolution=roi_lineshape_deconvolution;%
            end%
            
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
                
                if ~exist('supp_ti','var')
                    supp_ti=[num2str(data_set.acquno) '/' num2str(data_set.procno)];
                end
                figure(10);loglog(sc_pow10,val_pow10,'DisplayName',['S(+) ' supp_ti]);hold on  %    ));
                if  loop_over_spectra==20
                    figure(10);loglog(sc_pow10,val_pow10*sqrt(2),'k:','DisplayName',['1.41xS(+) NS=1']);hold on  %    ));
                    figure(10);loglog(sc_pow10,val_pow10*(2),'k--','DisplayName',['2xS(+) NS=1']);hold on  %    ));
                    
                end
                
                to_pt=round(min([size(sc_pow10_window,2) size(sc_pow10,2)])/2);
                v1=log(val_pow10_window(1,1:to_pt));
                v2=log(val_pow10(1:to_pt,1))';
                [~, pos_norm]=min((v2-v1));
                %
                
                drawnow;
                legend('Location','SouthWest');
                legend('off');
               
                if last%plot final figure
                    set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
                    y_dims=ylim;  ylim([1 y_dims(1,2)]);%set lower bound avoi showing long tail with very small values)
                    fig= figure(10);
                    set(fig,'PaperUnits' , 'centimeters');
                    %size of figure
                    width_fig=18;
                    width_fig=9;
                    width_fig=12;
                    height_fig=width_fig*2/3;
                     set(fig,'PaperOrientation','portrait');% portrait should be default...
                    set(fig,'PaperType', 'a4');% portrait should be default...
                    legend()
                    if ~exist('dataset_name','var')
                        dataset_name=['Dataset name : ' dataset];
                    end
                    title([dataset_name],'interpreter','none')
                    drawnow;
                    print('-dpdf',[ '.' filesep 'Results' filesep 'Final_comparison_snr_' data_set.expname '_' num2str(loop_over_spectra) '.pdf']);%here
                    figure(10);clf
                end
           
                    
            

        end
    end
end
%   pause

%print('-depsc','-tiff','-r600',['./Comparison_snr_' dataset '.eps']);%here
end


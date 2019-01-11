function ret=SIMPLE_DEMO()
%% function simple demo
path_of_experiments='./demo_nmr_data';
dataset='dj-caryophyllene_oxide';
exp_no=10;
exp_procno=1;

path_acqu =[ path_of_experiments filesep dataset filesep];

%[path_acqu]  = uigetdir(path_acqu);%remove this line to skip the user unterface to selec the spectrum
path_acqu=[path_acqu filesep];
full_path=[path_acqu  num2str(exp_no) filesep 'pdata' filesep num2str(exp_procno) filesep];
if exist(full_path ,'dir')
    
    disp(['Reading    spectrum  ' dataset  filesep num2str(exp_no) filesep 'pdata' filesep num2str(exp_procno)  ' ' ])
    data_set=read_data_bruker(path_acqu,exp_no,exp_procno);%read Bruker format
    
    disp(['Workin on spectrum  ' dataset filesep num2str(exp_no) filesep 'pdata' filesep num2str(exp_procno)  ' ' data_set.pulprog])
    
    mkdir('Results_paper_git')% in case does not exist
    
    %% set options values
    opt.fix_offset=1;
    opt.plot_results=1;%
    opt.fig_number=1;
    %opt.up_to_this_number_of_time_noise_level=5;
    
    %% determine noise level
    [data_set.noise_level, data_set.list_peaks, data_set.I0_offset, data_set.noise_levela , data_set.noise_leveln , data_set.noise_levelan, ...
        how_much_higher_than_noise_are_signals, where_determine_noise_level, sc_pow10, val_pow10, data_set.signal_shape] ...
        = get_noise_level_simple(data_set,opt);
    
    %% correct the spectrum
    data_set.spectrum= data_set.spectrum-data_set.I0_offset;
    
else
    error(['folder  ' full_path '  does not exist !' ])
end
ret=1;
end

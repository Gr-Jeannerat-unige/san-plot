function SAN_plot_from_bruker_data()
% SAN_plot_from_bruker_data()

upper_limit_of_expno_or_proco_to_reach=999;
msg=['Select the folder including the source for the SAN plot' 10 ];
msg=[msg 'For a single spectrum, point to .../expname/exp_number/pdata/proc_number ' 10 ];%type 1
msg=[msg 'For a serie of spectra with a given exp_number and processing numbers 1-' num2str(upper_limit_of_expno_or_proco_to_reach) ', point to : .../expname/exp_number ' 10 ];%type 2
msg=[msg 'For a serie of spectra with the exp_number 1-' num2str(upper_limit_of_expno_or_proco_to_reach) ' and processing 1, point to : .../expname ' 10 ];%type 3
width=700;
d = dialog('Position',[300 300 width 250],'Name','Generate SAN plot from Bruker spectra');

txt = uicontrol('Parent',d,'Style','text',...
    'Position',[20 80 width-(20+20) 80], 'String',msg);
w=70;

path_of_experiments=['.' filesep 'demo_nmr_data'];
dataset='';

base_path =[ path_of_experiments filesep dataset filesep];
counter=1;

while (counter==1) %loop until get data to generate plot
    
    [base_path]  = uigetdir(base_path);%('*','Point to fid or ser file');
    pos_find_sep=strfind(base_path,filesep);%filesep is the char of the folder separator of the sytem ('/'...)
    
    if ~(length(base_path)==pos_find_sep(1,size(pos_find_sep,2)))
        base_path=[base_path filesep];
        pos_find_sep=strfind(base_path,filesep);%filesep is the char of the folder separator of the sytem ('/'...)
    end
    %disp(['init : ' path_acqu2])
    m=size(pos_find_sep,2);
    
    delete(d)
    
    pos_find_pdata=strfind(base_path,'pdata');
    last_field_number=str2num(base_path(pos_find_sep(m-1)+1:pos_find_sep(m)-1));%try to get number
    
    %determine the type of folder (if single exp or series of procs or acqus)
    if size(pos_find_pdata,2)==0% no pdata in path
        if size(last_field_number,2)>0%last field is a number
            type=1;%multiprocs
        else
            type=3;%multiacqus
        end
    else
        %determine how many additional filesep after the pdata part
        n=size(find(pos_find_sep>pos_find_pdata),2);
        if n==2
            type=2;%single exp
        else
            type=1;%multiprocs
            %trim "pdata"
            base_path=base_path(1:pos_find_pdata-1);
            pos_find_sep=strfind(base_path,filesep);%filesep is the char of the folder separator of the sytem ('/'...)
            
            m=size(pos_find_sep,2);
            
            last_field_number=str2num(base_path(pos_find_sep(m-1)+1:pos_find_sep(m)-1));;%try to get number
            
        end
    end
    %disp(['inter: ' path_acqu2])
    
    m=size(pos_find_sep,2);
    if type==3
        list_acqus=[1:upper_limit_of_expno_or_proco_to_reach];
        list_procs=[1];
        base_path=base_path(1:pos_find_sep(m-0));
        
    end
    if type==1
        list_acqus =[last_field_number];
        list_procs=[1:upper_limit_of_expno_or_proco_to_reach];
        base_path=base_path(1:pos_find_sep(m-1));
        
    end
    if type==2
        wwh=2;
        list_acqus =[str2num(base_path(pos_find_sep(m-wwh-1)+1:pos_find_sep(m-wwh)-1))];
        list_procs=[last_field_number];
        base_path=base_path(1:pos_find_sep(m-3));
        
    end
    disp(['Base of the experimental data: ' base_path])
    
    for exp_no=list_acqus
        for exp_procno=list_procs
            full_path=[base_path  num2str(exp_no) filesep 'pdata' filesep num2str(exp_procno) filesep];
            if exist(full_path ,'dir')
                OK=1;
                
                disp(['Reading    spectrum  ' full_path  ' ' ])
                data_set=read_data_bruker(base_path,exp_no,exp_procno);%read Bruker format
                
                disp(['Workin on spectrum  ' full_path  ' ' data_set.pulprog])
                
                mkdir('Results_folder')% in case does not exist
                
                %% set options values
                opt.fix_offset=1;
                opt.plot_results=1;%
                opt.fig_number=100+counter;
                opt.up_to_this_number_of_time_noise_level=5;
                
                %% determine noise level
                [data_set.noise_level, data_set.list_peaks, data_set.I0_offset, data_set.noise_levela , data_set.noise_leveln , data_set.noise_levelan, ...
                    how_much_higher_than_noise_are_signals, where_determine_noise_level, sc_pow10, val_pow10, data_set.signal_shape] ...
                    = get_noise_level_simple(data_set,opt);
                
                %% correct the spectrum
                %
                %   for x = 1:10
                %       disp(x)
                %   end
                %
                data_set.spectrum= data_set.spectrum-data_set.I0_offset;
                
                counter=counter+1;
            end
        end
    end
    
end
end

function [ where_cut_stat ] = determine_position_measure_SNR(work_sp,opt,plot_results)
%determine the optimal position in the distribution to measure the noise
% i.e. where it varies the least about the central value about these
% different location tested

consider_list=[0.5:0.1:0.9];%testing beyond 0.9 requires to change the algorithm below

if nargin<3
    plot_results=0;
end
%UNTITLED2 determination position cut for SNR calculation
%   Detailed explanation goes here
fig_number_main=20000;
if isfield(opt,'fig_number')
    fig_number_main=20000+opt.fig_number;
end

test_list_of_values_where_cut_stat=0.01:0.01:0.99;
cutoff=zeros(2,size(test_list_of_values_where_cut_stat,2));
counter=1;
for cur_value=test_list_of_values_where_cut_stat
    val=(work_sp(round(size(work_sp,1)*cur_value),1))/(-simple_norminv((cur_value)/2));%level of signal at half the distribution of pos signals
    cutoff(1,counter)=cur_value;
    cutoff(2,counter)=val;
    counter=counter+1;
end
if plot_results
    figure(fig_number_main);clf;        plot(cutoff(1,:),cutoff(2,:),'r-'); hold on
end
min_ki=1e50;
for calc_this=consider_list
    delta=(consider_list(1,2)-consider_list(1,1))/2;
    from=calc_this-delta;
    to=calc_this+delta;
    list=find((cutoff(1,:)>from)&(cutoff(1,:)<to));%select position to take into account for averaging
    average_value=sum(abs(cutoff(2,list(1,:))/size(list,2)));%average
    if plot_results
        plot([cutoff(1,list(1,1)) cutoff(1,list(1,end))],average_value+[0 0],'k-')
    end
    ki_sq=sum(power(cutoff(2,list(1,:))-average_value,2));
    if plot_results
        plot(calc_this,average_value-sqrt(ki_sq),'rx');
        plot(calc_this,average_value+sqrt(ki_sq),'rx');
    end
    if ki_sq<min_ki
        min_ki=ki_sq;where_cut_stat=calc_this;level=average_value;
    end
end
if plot_results
    plot(where_cut_stat,level,'go');
    plot([0 1],level+[0 0],'g:');
    axis([0 1 0 average_value*2])
    title(['NSR at varions percentage of point distribution. Best for ' num2str(where_cut_stat)])
end
end


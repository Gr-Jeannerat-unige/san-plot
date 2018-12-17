function [sc_pow10,val_pow10]=interp_log_distrib(distr_wind,nb_pt_per_log10)
sca=1:size(distr_wind,2);
tot_pt_on_x_axis=round(nb_pt_per_log10*log10(size(distr_wind,2)));%number of points on the axis
end_of_log_scal=log10(sca(end));
start_of_log_scal=log10(sca(1));
sc_pow10=round(power(10,start_of_log_scal:((end_of_log_scal-start_of_log_scal)/(tot_pt_on_x_axis-1)):end_of_log_scal));
val_pow10=distr_wind(1,sc_pow10);
end


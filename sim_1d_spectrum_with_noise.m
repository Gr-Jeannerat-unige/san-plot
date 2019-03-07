function data = sim_1d_spectrum_with_noise(noise_level, lb, absolutevalue_mode,noise_in_fid,triplet)
plot_spectrum=1;
if nargin<1
    noise_level=1;% default noise level
end
if nargin<2
    lb=0.3;%default line width
end
if nargin<3
    absolutevalue_mode=0;%by default real spectrum (with 1: mc mode = magnitude or absolute value mode)
end
if nargin<4
    noise_in_fid=0;% default add noise in spectrum (better control of level... but does not passs the window function.)
end
if nargin<5
    triplet=0;% default add noise in spectrum (better control of level... but does not passs the window function.)
end
data=struct;

%% acquisition parameters
data.sfo1=500;%Larmor frequency
data.swp2=10;% spectral window in ppm
data.op1=5;% center of the windows
aq=5;%acquisition time in seconds

%% procesing parameters
nb_zero_filling=15;%(1, 3, 7, 15...)

swp=data.swp2;
sw=data.sfo1*data.swp2;
dw=1/sw;
th=0:dw/2:aq-dw/2;
t=th(1:2:end);

%% set signals's frequencies and intensities
nu=[];
height=[];
% add solvent
nu=[7.27*data.sfo1  ];
height=[0.2 ];
% add a triplet
nu=[nu 2.5*data.sfo1+[-7 0 7] ];
height=[height 1 2 1 ];

if ~triplet
height=[0 0 2 0 ];
end
% height=[1 0 0 0 ];
height=100*height;%the large signals will have a SINO of 100
%
% nu=[2.5*data.sfo1+[ 0 ] ];
% height=[200 ];

%% start acquisition
fid_full        =zeros(1,size(th,2));%(0)*th;
%generation of the fid / loop over signals
for loo_spi=1:size(nu,2)
    fid_full        =fid_full        +height(1,loo_spi)*exp(1i*2*pi*((th)*        nu(1,loo_spi))).*exp(-th*2*pi*lb);
end
detector_phase=0;
%apply phase of detector
fid_full        =(fid_full        *exp(1i*detector_phase/180*pi));

fid_sim=fid_full(1:2:end);% 1, 3, 5...
%%%fid_sim=awgn_dj(fid_sim,0);%add white gaussian noise 0 db
%  noise_db=-10*log(noise_level);
if noise_in_fid
    noise_db=-10*log(noise_level*280)/log(10);
    fid_sim=awgn_dj(real(fid_sim),noise_db)+1i*awgn_dj(imag(fid_sim),noise_db);%add white gaussian noise
else
    noise_db=-10*log(noise_level)/log(10);
    
end

fid_sim(1)=fid_sim(1)/2;
% apply zero filling and detection
fid_sim=[fid_sim zeros(1,nb_zero_filling*size(fid_sim,2))];

%apply Fourier transform
spectrum0=fftshift(fft((fid_sim)));%*exp(-1i*detector_phase/180*pi);
%spectrum0=fftshift(fft((fid_sim)))*exp(-1i*detector_phase/180*pi);
%spectrum0=fftshift(fft(([fid_sim zeros(1,nb_zero_filling*size(fid_sim))])))*exp(-1i*detector_phase/180*pi);
% rescale
spectrum0=spectrum0*max(max(height))/max(max(spectrum0));% rescale / could be done without need ot height... more elegants and more rigorous
% add noise 0 db : noise=1

disp(['Noise = ' num2str(noise_level) ' corresponds to -10log(Noise) = ' num2str(noise_db) ' db '])
if ~noise_in_fid
    spectrum0=awgn_dj(real(spectrum0),noise_db)+1i*awgn_dj(imag(spectrum0),noise_db);%add white gaussian noise 0 db
end
if absolutevalue_mode
    if absolutevalue_mode==2 %power mode
        spectrum0=      (real(spectrum0).^2+imag(spectrum0).^2);
    else %MC mode
        spectrum0=  sqrt(real(spectrum0).^2+imag(spectrum0).^2);
    end
else %pk mode (standard)
    spectrum0=real(spectrum0);
end

%% set scale
% set scale increment
incsw=swp/size(spectrum0,2);
% set scale
scale=data.op1+(-swp/2+incsw/2:incsw :swp/2-incsw/2);

plot_line_broadening=1.5;

if plot_spectrum
    
    figure(113);clf
    plot(scale,(spectrum0),'k-','Linewidth',plot_line_broadening);
end
data.td2=size(fid_sim,2);
data.tdeff2=size(fid_sim,2);
data.si2=size(spectrum0,2);
data.spectrum=spectrum0';
data.sw2h=sw;
data.wdw=1;
data.lb=lb;
data.expname='synthetic_spectrum';
data.acquno=1;
data.procno=1;
data.pulprog='zg';
data.ph_mod=1+absolutevalue_mode;

% %
% %             cur_f.Units='centimeters';
% %             tm=mod(inc_fig-10,4);
% %             cur_f.Position=[1+tm*16 1+10*(inc_fig-10-tm)/4 16 10];
% %
% %
% %             orient landscape
% %             print('-depsc','-tiff','-r600',[ 'Fig_demo_quad_1D_det_phase_' addtxt num2str(detector_phase) ' .eps']);%here
% %             if first_tp_non_zero
% %                 print('-depsc','-tiff','-r600',[ 'Fig_demo_quad_1D_det_phase_' addtxt num2str(detector_phase) '_first_p1_non_zero' num2str(first_tp_non_zero) '.eps']);%here
% %             end
% %             inc_fig=inc_fig+1;
end


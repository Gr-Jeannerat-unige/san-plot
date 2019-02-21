%% setting and reading spectrum 
file_id=fopen('/Volumes/lacie_case/nmr_data/nmrge500_3/data/winssinger/data/nmr/nmr/DC15-40/10/pdata/1/1r');
spectrum=fread(file_id,'int');
%% plotting spectrum for test and reference
figure(1);clf;plot(spectrum)
%% listing pos/neg points
index_of_positive_points=spectrum>0;
index_of_negative_points=spectrum<0;
%% generating the list of pos/neg points
Splus=abs(spectrum(index_of_positive_points));
Sminus=abs(spectrum(index_of_negative_points));
%% sorting points
Splus=sort(Splus,'descend');
Sminus=sort(Sminus,'descend');
%% plotting SAN plots of Splus and Sminus
figure(2);clf;
loglog(Splus);hold on
%loglog(Sminus)

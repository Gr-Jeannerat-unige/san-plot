%% setting and reading spectrum 
file_id=fopen('./dj-caryophyllene_oxide/12/pdata/4/1r');
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
Smax=Splus(1,1);
SigmaN=Splus(round(size(Splus,1)/2),1)/0.6744897495;
title(['Smax=' num2str(Smax) ', SigmaN=' num2str(SigmaN), ', SINO=' num2str(Smax/(2*SigmaN))])

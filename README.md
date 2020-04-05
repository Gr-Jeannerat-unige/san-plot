# san-plot
Software for Signal-Artifacts-Noise (SAN) plot.

This repository contains:

1) The program used to generate the figures of the [paper](https://onlinelibrary.wiley.com/doi/full/10.1002/mrc.4882) : 

program_used_to_prepare_san_paper.m (probably only works with matlab on mac)

2) The folder including the data used for the demo/article : demo_nmr_data

3) The very simple demo program : SIMPLE_DEMO.m (working using matlab and octave 4.4)

4) The function to generate SAN plots from Bruker data : SAN_plot_from_bruker_data.m

Note that the interface to select the folder may not work with some versions of Octave.

When starting, the programs directly opens a path-selection window.

Select the folder including the data you want to use to plot the SAN plot.

For a single spectrum, point to:

.../expname/exp_number/pdata/proc_number 

For a serie of spectra with a given exp_number and processing numbers 1-999, point to: 

.../expname/exp_number

For a serie of spectra with the exp_number 1-999  and processing 1, point to:

 .../expname 


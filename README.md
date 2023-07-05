# 3C273DiskWindPaper
Code and related files for Long+2023

The model (lightly documented in the code, more robustly documented in the paper) is contained as the Julia module [`DiskWind`](DiskWind). A sample Python fitting script that uses the model is contained in [`PTemceeFit.py`](PTemceeFit.py) (which is run on HPC resources with [`run.py`](run.py)). The code uses a convention for position angle that is 90 degrees offset from the standard astronomy definition, such that the values obtained from fitting are 90 degrees higher than they should be, and fits the mass in terms of $3\times 10^8 M\odot$ to compare with the rough value previously determined for 3C 273. All of our fit results (chunked for job size purposes) using the model and this script are available as pickle files (`jPyPTEmceeVar*.p`).

Scripts that reproduce the figures are available in [`figures.py`](figures.py) and [`figures.jl`](figures.jl). To create the reverberation mapping figure some additional processing must be done with [`RMDelaySummit.jl`](RMDelaySummit.jl) which outputs products `RMVars_*.jld` depending upon choice of parameters. 

The GRAVITY data used to fit the model are contained in the Pickle file [`3c273_juljanmarmay_append_gilles_specirf_wide_v6.p`](3c273_juljanmarmay_append_gilles_specirf_wide_v6.p) and were obtained as described in [GRAVITY+ (2018)](https://www.nature.com/articles/s41586-018-0731-9): "The data were obtained at the VLTI of the European Southern Observatory (ESO), Paranal, Chile, and are available on the ESO archive (http://archive.eso.org/eso/eso_archive_main.html) under programme IDs 099.B-0606, 0100.B-0582 and 0101.B-0255." 

If you have any questions please reach out as I am more than happy to help/explain more, and if you find this work useful please cite our paper.

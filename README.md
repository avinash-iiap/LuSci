# Lunar Scintillometer
This code repository is the core software pipeline for analyzing the RAW data coming from a six-channel lunar scintillometer as described in the paper, 'Development of a Lunar Scintillometer as part of the National Large Optical Telescope Site Survey' which is submitted to the journal, 'Experimental Astronomy'. The pipeline has the LabView code for acquiring data from the instrument, and an automated bash script to analyze the RAW voltage data and generate the atmospheric seeing values. The code is divided into two parts:

1. A **stand-alone code** which processes a single LVM file from a single night and creates the vertical seeing distribution for the same.

2. A **campaign mode** which processes all the LVM files in a specified input folder, creates the vertical seeing distribution for all the nights and generates an overall analysis for the whole campaign.

## Stand-alone mode
**lusci_run1.sh** is the bash script which runs the stand-alone pipeline which processes the data for a single night. It runs the following codes one after the other:

1. **Main_code/voltage_to_covariance_500samp_v1.m** is a signal pre-processing package which takes the LabView LVM file as the input file, generates the covariances between the different baselines and writes it to a file named **lusci1_yyyy-dd-mm.cov** (yyyy-dd-mm is replaced by the night at which the data was taken). The .cov file acts as the input file for the profrest.pro IDL package.

2. **profrest.pro (NOT included in the github package)** is the covariance to optical turbulence profile (OTP) converter. It takes the lusci1_yyyy-dd-mm.cov as the input file, takes a set of input parameters from lusci_hanle.par (Hanle, Leh is the site at which the instrument was tested with a campaign run for 2 years), and generates a **lusci1_yyyy-dd-mm.tp** file which has the atmospheric structure constants and the turbulence integrals for all the pivot point altitudes. The code is available as a tarball with dependencies by the name, **code4.tar.gz** from the original author, Prof. A Tokovinin of CTIO. **IMPORTANT: The code4 folder should be present in the /Main_code subfolder for the stand-alone code to work. Prof. A. Tokovinin should be contacted directly to get the relevant code**. For more details on the restoration algorithm and how profrest.pro works, kindly check the following links:
 - http://www.ctio.noao.edu/~atokovin/profiler/code4.pdf
 - http://www.ctio.noao.edu/~atokovin/profiler/restor3.pdf
  
3. **Main_code/report_compiled_v1.m** is a plotting package which displays the turbulence statistics for each night. It relies on the **lusci1_yyyy-dd-mm.tp** files created by **profrest.pro** and the .mat file with the RAW voltage parameters created by **Main_code/voltage_to_covariance_500samp_v1.m** to generate the plots. For more details on the plots created, please refer to the comments in the code header.

### Input parameter files and Dependencies

1. **lusci_params1.txt** contains the input parameters used to point to the file paths and generate the covariance file.
 - input_file	'/home/cyrano/LuSci/Hanle/obs_05042015.lvm'
 - output_folder	'/home/cyrano/LuSci/cov_outputs/automated_new/2015_04_Apr/'
 - File index	1
 - IDL_parameters	'/home/cyrano/LuSci/lusci_hanle.par'
  Dark		'/home/cyrano/LuSci/Hanle/dark_05042015.lvm'
Sky		''
Samples		500

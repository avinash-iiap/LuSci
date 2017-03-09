# Lunar Scintillometer
This code repository is the core software pipeline for analyzing the RAW data coming from a six-channel lunar scintillometer as described in the paper, 'Development of a Lunar Scintillometer as part of the National Large Optical Telescope Site Survey' which is submitted to the journal, 'Experimental Astronomy'. The pipeline has the LabView code for acquiring data from the instrument, and an automated bash script to analyze the RAW voltage data and generate the atmospheric seeing values. 

**LuSci DAQ.vi** is the LabView package used to acquire data from the NI USB-6210 DAQ. Make sure the DAQ device is connected to the USB port before opening this file. It consists of a simple six-channel acquisition block set to run at 500 Hz and a power spectrum analysis block to see whether any EMI interference is present in the incoming data (in real-time). 

The software pipeline is divided into two parts:

1. A **stand-alone code** which processes a single LVM file from a single night and creates the vertical seeing distribution for the same.

2. A **campaign mode** which processes all the LVM files in a specified input folder, creates the vertical seeing distribution for all the nights and generates an overall analysis for the whole campaign.

## Stand-alone mode
**lusci_run1.sh** is the bash script which runs the stand-alone pipeline which processes the data for a single night. It runs the following codes one after the other:

1. **Main_code/voltage_to_covariance_500samp_v1.m** is a signal pre-processing package which takes the LabView LVM file as the input file, generates the covariances between the different baselines and writes it to a file named **lusci1_yyyy-dd-mm.cov** (yyyy-dd-mm is replaced by the night at which the data was taken). The .cov file acts as the input file for the profrest.pro IDL package.

2. **profrest.pro (NOT included in the github package)** is the covariance to optical turbulence profile (OTP) converter. It takes the lusci1_yyyy-dd-mm.cov as the input file, takes a set of input parameters from lusci_hanle.par (Hanle, Leh is the site at which the instrument was tested with a campaign run for 2 years), and generates a **lusci1_yyyy-dd-mm.tp** file which has the atmospheric structure constants and the turbulence integrals for all the pivot point altitudes. The code is available as a tarball with dependencies by the name, **code4.tar.gz** from the original author, Prof. A Tokovinin of CTIO (email:atokovinin@ctio.noao.edu). **IMPORTANT: The code4 folder should be present in the /Main_code subfolder for the stand-alone code to work. Prof. A. Tokovinin should be contacted directly to get the relevant code**. For more details on the restoration algorithm and how profrest.pro works, kindly check the following links:
 - http://www.ctio.noao.edu/~atokovin/profiler/code4.pdf
 - http://www.ctio.noao.edu/~atokovin/profiler/restor3.pdf
  
3. **Main_code/report_compiled_v1.m** is a plotting package which displays the turbulence statistics for each night. It relies on the **lusci1_yyyy-dd-mm.tp** files created by **profrest.pro** and the .mat file with the RAW voltage parameters created by **Main_code/voltage_to_covariance_500samp_v1.m** to generate the plots. For more details on the plots created, please refer to the comments in the code header.

### Input parameter files

1. **lusci_params1.txt** contains the input parameters used to point to the file paths and generate the covariance file.
 1. input_file	- Full path to the input LabView LVM file which has the RAW voltages from the six channels. The LVM file should be as specified under LabView software standards.
 2. output_folder	- The folder to which the .cov,.tp and the results file are to be exported.
 3. File index	- Default value is 1. If you have more than one file per night, increment this value for every file to ensure that the output files within each night don't get overwritten.
 4. IDL_parameters	- The parameter file used by the IDL package to compute OTP.
 5. Dark	- Location of the dark file. The dark file should the RAW voltages from the six photodiodes taken through 'LuSci DAQ.vi' while the photodiodes are completely covered from light. It should be of at least 5-10 minutes duration to ensure that the filter capacitors are fully discharged. If not used, use '' here.
 6. Sky	- Location of the sky file. The sky file should the RAW voltages from the six photodiodes taken through 'LuSci DAQ.vi' while the photodiodes are exposed to an empty patch of the sky away from the moon. It should be of at least 5-10 minutes duration to ensure that the filter capacitors are fully discharged, and to ensure that no anomalies in that patch of the sky contribute to incorrect data. If not used, use '' here.
 7. Samples	- The sampling rate (in samples/second) used for acquiring the RAW voltage data.
 
2. **lusci_hanle.par** contains the input parameters for the IDL package including the geographical location and pivot point altitudes. Most of the parameters are commented and are self explanatory. For more details, check the IDL package manuals mentioned under the profrest.pro section. You can rename this file, but make sure to change it in **lusci_params1.txt** and **lusci_params_complete.txt**.

## Campaign mode
The campaign mode wraps around the stand-alone mode. The campaign mode accepts the input folder where all the LVM files corresponding to the raw voltage data from different nights are located, and executes the stand-alone analysis of each night after setting the input parameters for each night of observation. **lusci_run_complete.sh** is the bash script which runs the complete campaign pipeline which processes the data for the whole campaign. 

### Format of the input LVM files for campaign mode
- Observation files - obs_ddmmyyyy.lvm
- Dark files - dark_ddmmyyyy.lvm
- Sky background files - sky_ddmmyyyy.lvm
where the 'ddmmyyyy' is replaced by the observation date. The LVM files *of each night* should be compressed into a single zip archive and all the zip archives of the different nights (for which the observation was taken) placed in the input folder (one of the parameters to **lusci_params_complete.txt**). In a single zip archive, there should compulsarily be an observation file but the dark and sky files are optional. There is no naming convention for the zip files but there should be no other files in the input folder. This format is not applicable for stand-alone mode where the .lvm files should be used as is. For an example of how an input file should be for the campaign mode, download this [file](https://drive.google.com/open?id=0BzDZNA262Mq0Y2dyZ2NQTEU1Vmc). This is a sample observation file taken from Hanle for a few minutes on 18/09/2013.
**NOTE**: Multiple observation files on a single night are not allowed in campaign mode unlike the stand-alone mode (where multiple output files can be changed with the File_index parameter in **lusci_params1.txt**). Combine the LVM files in ascending order of their time of observation and create a new LVM file (an automated code for combining LVM files will be posted soon).

### Input parameter file
**lusci_params_complete.txt** contains the input parameters used to point to the folder where all input files are located and to the folder where the output results are to be generated.
1. input_folder	- the full path of the folder where the compressed zip archive(s) of the input LVM files are located.
2. output_folder	- the full path of the folder where the output .cov, .tp and the results for each night will be generated and stored in folders named after each month.
3. final_analysis	- the full path of the folder where the full campaign analysis results are stored. 
4. IDL_parameters	- Same as in lusci_params1.txt
5. Samples	- Same as in lusci_params1.txt

**lusci_run_complete.sh** clears the output folders, finds the number of zip files, extracts each zip file iteratively and runs the **lusci_run1.sh** for each night of observation. It deletes each LVM file after processing is done for that night. After all the files of the campaign are processed, **Main_code/compiled_seeing.m** is invoked to generate the final campaign results in the folder specified by 'final_analysis' parameter in **lusci_params_complete.txt**. See the header comments in **compiled_seeing.m** to understand the nature of the results.

The code is licensed under The GNU General Public License v3.0, and can be freely modified and distributed. I acknowledge [M. A. Hopcroft](https://in.mathworks.com/matlabcentral/profile/authors/633845-m-a-hopcroft) for the MATLAB package which imports LabView standard LVM files.

Author Name : Avinash Surendran
Affiliation : Indian Institute of Astrophysics, Bangalore, India
Email : asurendran@iiap.res.in

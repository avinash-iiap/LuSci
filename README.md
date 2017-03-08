# Lunar Scintillometer
This code repository is the core software pipeline for analyzing the RAW data coming from a six-channel lunar scintillometer as described in the paper, 'Development of a Lunar Scintillometer as part of the National Large Optical Telescope Site Survey' which is submitted to the journal, 'Experimental Astronomy'. The pipeline has the LabView code for acquiring data from the instrument, and an automated bash script to analyze the RAW voltage data and generate the atmospheric seeing values. The code is divided into two parts:

1. A **stand-alone code** which processes a single LVM file from a single night and creates the vertical seeing distribution for the same.

2. A **campaign mode** which processes all the LVM files in a specified input folder, creates the vertical seeing distribution for all the nights and generates an overall analysis for the whole campaign.

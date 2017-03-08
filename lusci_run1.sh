#!/bin/bash

# LuSci single-night mode automation bash script. 
# Executes 'Voltage to Covariance' MATLAB script, reads the parameters written into TMP.txt by the script and passes them as arguments on to the IDL program. The single-night plotting script is executed after the IDL prohram computes the OTP (Optical Turbulence Profile).

cp lusci_params1.txt ${PWD}/Main_code/lusci_params1.txt
cd Main_code
sudo matlab -nodesktop -nosplash -r "voltage_to_covariance_500samp_v1;quit;"
rm -f lusci_params1.txt
read -a d -d '\n' < TMP.txt
# Following line selects the minimum size of cov file to continue processing
if [ ${d[3]} -le 10 ]
then
echo $"-----------------Usable data is too small, exiting the program-------------------"
exit
fi
cd ${PWD}/code4
sudo idl -e "profrest" -args ${d[0]} ${d[1]} ${d[2]}
cd ..
sudo matlab -nodesktop -nosplash -r "report_compiled_v1;quit;"
rm -f TMP.txt




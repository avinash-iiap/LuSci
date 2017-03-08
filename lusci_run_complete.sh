#!/bin/bash

# LuSci campaign mode automation bash script.
# Reads the lusci_params_complete.txt for the input folder where the .lvm file resides, and output folders where the intra-night results and the complete campaign analysis is to be stored. It wraps over the lusci_run1.sh which is the single-night mode automation script. For the correct working of this script, there should only be the observation, dark and sky files of the relevant nights in the input folder. All lvm files need to be archived as a single zip file or in multiple zip files. The files inside the archive should be named as follows:
# Observation files - obs_ddmmyyyy.lvm
# Dark files - dark_ddmmyyyy.lvm
# Sky background files - sky_ddmmyyyy.lvm
# where ddmmyyyy should be replaced by the date of the observation in that specific format

readarray params < 'lusci_params_complete.txt'
home=${PWD}
input=${params[0]#*\'}
input=${input%\'*}
output=${params[1]#*\'}
output=${output%\'*}
finale=${params[2]#*\'}
finale=${finale%\'*}
cd "${output}"
find . -name "*" -exec rm -rf {} \;
echo "Old results in the output folder deleted"
cd "${finale}"
find . -name "*" -exec rm -rf {} \;
echo "Old final analysis results deleted"
#cd "${home}/Main_code"
#find . -name "*.log" -exec rm -rf {} \;
cd "${input}"
find . -name "*.lvm" -exec rm -rf {} \;
echo "Old lvm files deleted"
if [ -e "${home}/Main_code/compiled_seeing.csv" ]; then
	rm -rf "${home}/Main_code/compiled_seeing.txt" 
	rm -rf "${home}/Main_code/compiled_seeing.csv"
	echo "Old compiled seeing files found and deleted"
fi
find . -type f -name '*.zip' -print0 | while IFS= read -r -d '' file; do
    (	cd "${input}"
	printf '%s\n' "Processing file $file"
	unzip -j "$file"
	file_check=$(ls obs*.lvm 2> /dev/null | wc -l)
	if [ $file_check -gt 0 ]; then
		pattern="obs*.lvm"
		file=( $pattern )
		printf "%b" "input_file\t'${PWD}/${file}'\n" > ${home}/lusci_params1.txt	
	else
		echo "Observation file not found; Skipping the current night"
	fi
	datestr=${file:4:8};
	fldrname=$(date -d ${datestr:4:4}${datestr:2:2}${datestr:0:2} +%Y_%m_%b)
	echo "Date - ${datestr}"
	printf "%b" "output_folder\t'${output}${fldrname}/'\n" >> ${home}/lusci_params1.txt
	printf "%b" "File index\t1\n" >> ${home}/lusci_params1.txt
	printf "%b" "${params[3]}" >> ${home}/lusci_params1.txt
	file_check=$(ls dark*.lvm 2> /dev/null | wc -l)
	if [ $file_check -gt 0 ]; then
		pattern="dark*.lvm"
		file=( $pattern )
		printf "%b" "Dark\t\t'${PWD}/${file}'\n" >> ${home}/lusci_params1.txt
		echo "Dark file found"	
       	else
		echo "Dark file not found"
		printf "%b" "Dark\t\t''\n" >> ${home}/lusci_params1.txt
	fi

	file_check=$(ls sky*.lvm 2> /dev/null | wc -l)
	if [ $file_check -gt 0 ]; then
		pattern="sky*.lvm"
		file=( $pattern )
		printf "%b" "Sky\t\t'${PWD}/${file}'\n" >> ${home}/lusci_params1.txt	
		echo "Sky file found"
       	else
		echo "Sky file not found"
		printf "%b" "Sky\t\t''\n" >> ${home}/lusci_params1.txt
	fi
	printf "%b" "${params[4]}" >> ${home}/lusci_params1.txt
	file_check=$(ls obs*.lvm 2> /dev/null | wc -l)
	if [ $file_check -gt 0 ]; then
		cd "${home}"
		sudo ./lusci_run1.sh
	fi
	cd "${input}"
	find . -name "*.lvm" -exec rm -rf {} \; )
done
wait
if [ ! -d "$finale" ]; then
	mkdir "${finale}" 
fi
find "${output}" -wholename "*seeing_data.txt" -print0 | xargs -0 sudo cp -t "${finale}"
cd ${home}/Main_code
sed -e 's/\s\+/,/g' compiled_seeing.txt > compiled_seeing.csv
sudo cp "${home}/Main_code/compiled_seeing.m" "${finale}"
sudo cp "${home}/Main_code/compiled_seeing.csv" "${finale}"
cd "${finale}"
sudo matlab -nodesktop -nosplash -r "compiled_seeing;quit;"

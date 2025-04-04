#!/bin/bash
# Make sure the script is executable: chmod +x runJob.sh
cd /afs/cern.ch/user/m/mcrucian/private/displaced_muons/CMSSW_13_3_0/src/DisplacedMuons-FrameWork/Ntuplizer/test
cmsenv

# Get command line arguments: input directory and output file name
input_dir=$1
out_file=$2

# Define log file (this can be customized)
logfile="log_MC_${input_dir}_MiniAOD.log"

echo "Running cmsRun for input directory: ${input_dir}"
start_time=$(date +%s)

# Run cmsRun with your configuration file; adjust if needed
cmsRun MC_MiniAOD_runNtuplizer_cfg_test.py -input_dir ${input_dir} -out_file ${out_file} &> ${logfile}

end_time=$(date +%s)
echo "Time taken: $((end_time - start_time)) seconds"

# Copy the output file to the EOS output directory if it exists
output_dir="/eos/home-m/mcrucian/SWAN_projects/DisplacedMuons/ntuples"
if [ -f "${out_file}" ]; then
    cp ${out_file} ${output_dir}
    echo "${out_file} copied to ${output_dir}"
else
    echo "Error: ${out_file} not found."
fi

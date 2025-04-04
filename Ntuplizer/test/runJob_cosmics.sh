#!/bin/bash
# Make sure the script is executable: chmod +x runJob_cosmics.sh
cd /afs/cern.ch/user/m/mcrucian/private/displaced_muons/CMSSW_13_3_0/src/DisplacedMuons-FrameWork/Ntuplizer/test
cmsenv

# Define log file (this can be customized)
logfile="log_Cosmics_AOD.log"

echo "Running cmsRun for Cosmics..."
start_time=$(date +%s)

# Run cmsRun with your configuration file
cmsRun Cosmics_runNtuplizer_AOD_cfg.py &> ${logfile}

end_time=$(date +%s)
echo "Time taken: $((end_time - start_time)) seconds"

# Copy the output file to the EOS output directory if it exists
output_dir="/eos/home-m/mcrucian/SWAN_projects/DisplacedMuons/ntuples"
if [ -f "ntuples.root" ]; then
    cp ntuples.root ${output_dir}
    echo "ntuples.root copied to ${output_dir}"
else
    echo "Error: ntuples.root not found."
fi
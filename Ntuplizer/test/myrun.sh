#!/bin/bash

cmsenv

output_dir="/eos/home-m/mcrucian/SWAN_projects/DisplacedMuons/"

echo "Running cmsRun for HTo2ZdTo2mu2x_MZd-40_Epsilon-2e-08..."
cmsRun MC_MiniAOD_runNtuplizer_cfg.py -input_dir HTo2ZdTo2mu2x_MZd-40_Epsilon-2e-08 -out_file MC_Ntuples_MZd-40_Epsilon-2e-08.root &>log_MC_MZd-40_Epsilon-2e-08_MiniAOD.log

echo "Running cmsRun for HTo2ZdTo2mu2x_MZd-40_Epsilon-8e-08..."
cmsRun MC_MiniAOD_runNtuplizer_cfg.py -input_dir HTo2ZdTo2mu2x_MZd-40_Epsilon-8e-08 -out_file MC_Ntuples_MZd-40_Epsilon-8e-08.root &>log_MC_MZd-40_Epsilon-8e-08_MiniAOD.log

echo "Running cmsRun for HTo2ZdTo2mu2x_MZd-40_Epsilon-5e-09..."
cmsRun MC_MiniAOD_runNtuplizer_cfg.py -input_dir HTo2ZdTo2mu2x_MZd-40_Epsilon-5e-09 -out_file MC_Ntuples_MZd-40_Epsilon-5e-09.root &>log_MC_MZd-40_Epsilon-5e-09_MiniAOD.log

echo "Running cmsRun for Cosmics..."
cmsRun Cosmics_runNtuplizer_AOD_cfg.py &> log_Cosmics_AOD.log

echo "Checking and copying ntuples.root..."
if [ -f "ntuples.root" ]; then
    cp ntuples.root "$output_dir"
    echo "ntuples.root copied to $output_dir"
else
    echo "Error: ntuples.root not found."
fi

echo "Checking and copying MC_Ntuples_MZd-40_Epsilon-2e-08.root..."
if [ -f "MC_Ntuples_MZd-40_Epsilon-2e-08.root" ]; then
    cp MC_Ntuples_MZd-40_Epsilon-2e-08.root "$output_dir"
else
    echo "Error: MC_Ntuples_MZd-40_Epsilon-2e-08.root not found."
fi

echo "Checking and copying MC_Ntuples_MZd-40_Epsilon-8e-08.root..."
if [ -f "MC_Ntuples_MZd-40_Epsilon-8e-08.root" ]; then
    cp MC_Ntuples_MZd-40_Epsilon-8e-08.root "$output_dir"
else
    echo "Error: MC_Ntuples_MZd-40_Epsilon-8e-08.root not found."
fi

echo "Checking and copying MC_Ntuples_MZd-40_Epsilon-5e-09.root..."
if [ -f "MC_Ntuples_MZd-40_Epsilon-5e-09.root" ]; then
    cp MC_Ntuples_MZd-40_Epsilon-5e-09.root "$output_dir"
else
    echo "Error: MC_Ntuples_MZd-40_Epsilon-5e-09.root not found."
fi

echo "Script execution completed."
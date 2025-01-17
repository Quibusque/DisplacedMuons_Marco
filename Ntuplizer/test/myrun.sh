#!/bin/bash

cmsenv

output_dir="/eos/home-m/mcrucian/SWAN_projects/DisplacedMuons/"

cmsRun MC_MiniAOD_runNtuplizer_cfg.py -input_dir HTo2ZdTo2mu2x_MZd-40_Epsilon-2e-08 -out_file MC_Ntuples_MZd-40_Epsilon-2e-08.root &>log_MC_MZd-40_Epsilon-2e-08_MiniAOD.log

cmsRun MC_MiniAOD_runNtuplizer_cfg.py -input_dir HTo2ZdTo2mu2x_MZd-40_Epsilon-8e-08 -out_file MC_Ntuples_MZd-40_Epsilon-8e-08.root &>log_MC_MZd-40_Epsilon-8e-08_MiniAOD.log

cmsRun MC_MiniAOD_runNtuplizer_cfg.py -input_dir HTo2ZdTo2mu2x_MZd-40_Epsilon-5e-09 -out_file MC_Ntuples_MZd-40_Epsilon-5e-09.root &>log_MC_MZd-40_Epsilon-5e-09_MiniAOD.log

cmsRun Cosmics_runNtuplizer_AOD_cfg.py &> log_Cosmics_AOD.log

if [ -f "ntuples.root" ]; then
    cp ntuples.root "$output_dir"
else
    echo "Error: ntuples.root not found."
fi

if [ -f "MC_Ntuples_MZd-40_Epsilon-2e-08.root" ]; then
    cp MC_Ntuples_MZd-40_Epsilon-2e-08.root "$output_dir"
else
    echo "Error: MC_Ntuples_MZd-40_Epsilon-2e-08.root not found."
fi

if [ -f "MC_Ntuples_MZd-40_Epsilon-8e-08.root" ]; then
    cp MC_Ntuples_MZd-40_Epsilon-8e-08.root "$output_dir"
else
    echo "Error: MC_Ntuples_MZd-40_Epsilon-8e-08.root not found."
fi

if [ -f "MC_Ntuples_MZd-40_Epsilon-5e-09.root" ]; then
    cp MC_Ntuples_MZd-40_Epsilon-5e-09.root "$output_dir"
else
    echo "Error: MC_Ntuples_MZd-40_Epsilon-5e-09.root not found."
fi
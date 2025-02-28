# Displaced Muons Framework

Taken from [DisplacedMuons-FrameWork](https://github.com/24LopezR/DisplacedMuons-FrameWork)

## Instructions for Installing and Producing Cosmic Data/MC Plots

### Installing

I recommend using `CMSSW_13_2_0`, but any later release should also work. For installing the code, follow these instructions:

```bash
cmsrel CMSSW_13_2_0
cd CMSSW_13_2_0/src
cmsenv
git clone git@github.com:24LopezR/DisplacedMuons-FrameWork.git
scram b -j8
```

### Gen matching - How I Run

At the moment, I am focusing on running only `ntuplizer_test` which considers only DSA muons for MC signals of displaced muons (ntuplizer runs also for the cosmic rays).

Example:

```bash
cmsRun MC_MiniAOD_runNtuplizer_cfg_test.py -input_dir HTo2ZdTo2mu2x_MZd-40_Epsilon-5e-09 -out_file test.root &> myLog.log
```

In the file `MC_MiniAOD_runNtuplizer_cfg_test.py`, you can edit the variable `nEvents` to run over a different number of events. Use `-1` to run over all events.


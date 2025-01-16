# Displaced Muons Framework

Taken from https://github.com/24LopezR/DisplacedMuons-FrameWork

## Instructions for installing and producing Cosmic Data/MC plots
### Installing
I recommend using `CMSSW_13_2_0`, but any later release should also work.
For installing the code, follow these instructions:

    cmsrel CMSSW_13_2_0
    cd CMSSW_13_2_0/src
    cmsenv
    git clone git@github.com:24LopezR/DisplacedMuons-FrameWork.git
    scram b -j8
Note that for computing the efficiencies and producing the plots we will use only the "Ntuplizer" package. The other two ("Analyzer" and "MiniAOD_prod") can be ignored / deleted.
### Producing Ntuples
To produce the Ntuples, one can start from AOD or MiniAOD. Depending if you are running on AOD/MiniAOD and Data/MC, you have to change the configuration of the ntuplizer. The four combinations are in `DisplacedMuons-FrameWork/Ntuplizer/python`
Then to run the ntuplizer you have to `cmsRun` any of the scripts you find in `DisplacedMuons-FrameWork/Ntuplizer/test` .
An example is summarized in these instructions:

    cd DisplacedMuons-FrameWork/Ntuplizer/test
    
    # Check listOfFiles in 'CosmicsMC_MiniAOD_runNtuplizer_cfg.py':
    #     listOfFiles = ['']
    # Check configuration in 'CosmicsMC_MiniAOD_runNtuplizer_cfg.py':
    #     process.load("DisplacedMuons-FrameWork.Ntuplizer.CosmicsMC_ntuples_MiniAOD_cfi")
    
    cmsRun CosmicsMC_MiniAOD_runNtuplizer_cfg.py &> log_CosmicsMC_MiniAOD.log
### Plotting efficiencies
For all the plots I use the script `plot_efficiencies.py`.
Usage:

        usage: plot_efficiencies.py [-h]
                                [--var [{dmu_dsa_dxy,dmu_dsa_pt,dmu_dsa_eta,dmu_dgl_dxy,dmu_dgl_dz} [{dmu_dsa_dxy,dmu_dsa_pt,dmu_dsa_eta,dmu_dgl_dxy,dmu_dgl_dz} ...]]]
                                [--mcfile MCFILE] [--datafile DATAFILE]
    
    optional arguments:
      -h, --help            show this help message and exit
      --var [{dmu_dsa_dxy,dmu_dsa_pt,dmu_dsa_eta,dmu_dgl_dxy,dmu_dgl_dz} [{dmu_dsa_dxy,dmu_dsa_pt,dmu_dsa_eta,dmu_dgl_dxy,dmu_dgl_dz} ...]]
                            Variable(s) to plot
      --mcfile MCFILE       MC Ntuple file
      --datafile DATAFILE   Data Ntuple file
The options I used are:

    python3 plot_efficiencies.py --var dmu_dsa_dxy
    python3 plot_efficiencies.py --var dmu_dsa_eta
    python3 plot_efficiencies.py --var dmu_dsa_pt
    python3 plot_efficiencies.py --var dmu_dgl_dxy dmu_dgl_dz
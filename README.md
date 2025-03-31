# pac-plv-gmbi
Arun Asthagiri

Scripts for PLV and PAC analyses. 

PAC_scripts/
1. analysis_GAMpac_surrogate.m: comodulogram using Canolty PAC modulation index sweeping both phase frequencies and amplitude frequencies.
2. analysis_PAC_phase_bins.m: plotting distribution of gamma amplitude over theta phases. Fixed theta frequency (defined per song), multiple gamma frequencies.
3. analysis_window_PAC_fois.m: sliding window PAC for hypothesized theta and gamma frequencies. (Used for main gMBI PAC analysis). 

PLV_script/
1. analysis_PLV.m: phase-locking value over frequencies using Morlet wavelet decomposition.

stargate_scripts/ # used for gMBI analysis across participants, sessions, conditions, songs
1. params_gMBI.m: hyperparameters for included participants
2. runPACPipeline.m: call PAC_scripts files within this function to integrate with gMBI analysis.
3. runPLVPipeline.m: call PLV_script file within this function to integrate with gMBI ananalysis (effectively same as runPACPipeline.m). 

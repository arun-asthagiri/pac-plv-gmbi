% params_gMBI.m
% parameters for gamma MBI EEG analysis

 %{
    BMK - March 2024
    Based on original psd/topo analysis by Ed Large

"/work/mindlab/Programs/eeglab2021.1/plugins/Fieldtrip-lite20210601/"
from your path, see http://bit.ly/2SPPjUS
 %}

%% path parameters
params.homedir = '/home/b.kubit/gammaMBI/Ring_EEG/';
params.path2data = '/work/mindlab/Projects/GammaMBI/EEG/EEG_Data/';
params.eeglabdir = fullfile('/work/mindlab/Programs/eeglab2021.1'); %params.eeglabdir = fullfile(params.homedir,'/eeglab-eeglab2021.1');
addpath(genpath(params.eeglabdir));
rmpath(genpath(fullfile(params.eeglabdir,'plugins/Fieldtrip-lite20210601'))); %idk what this is,but causes problems with signal toolbox code
params.path2electLocation = strcat(params.eeglabdir,'/plugins/dipfit4.3/standard_BEM/elec/standard_1005.elc'); %strcat(params.eeglabdir,'/plugins/dipfit5.3/standard_BESA/standard-10-5-cap385.elp');  %);

%% pick channels for analysis (not preprocessing)
chans2use=[];%["O1","OZ","O2"]; %[] = all channelsl ["O1","OZ","O2","P7","P8"]

%% choose preprocessing scheme
preproc = ''; %options: '', 'p1', 'p2',
%^ files names of preprocessed song files
if strcmp(preproc,'')
    %Ed changes (nothing appended to file names)
    params.anals.winsecs = 120;
    [b,a] = nozFilter(1,1); %Filter for local gain calc, Nozaraden-like (1,3)
    params.anals.chans = [1:30];
    params.preproc.refType = 'tps'; %tps-TP9 + TP10; avg average ref; dryeeg A2
elseif strcmp(preproc,'p1')
    %orig parameters (EVT1 test)
    params.anals.winsecs = 30;
    [b,a] = nozFilter(1,2); %Filter for local gain calc, Nozaraden-like
    params.anals.chans = [1:30];
    params.preproc.refType = 'tps'; %tps-TP9 + TP10; avg average ref
elseif strcmp(preproc,'p2')
    %avg reference instead of TP9+10
    params.anals.winsecs = 120;
    [b,a] = nozFilter(1,3); %Filter for local gain calc, Nozaraden-like
    params.anals.chans = [1:32];
    params.preproc.refType = 'avg'; %tps-TP9 + TP10; avg average ref
elseif strcmp(preproc,'p4')
    %Ed changes trying diff filter for each freq. same preproc as ''
    params.anals.winsecs = 120; %used for theta and delta
    [b,a] = nozFilter(1,3); %Filter for local gain calc, Nozaraden-like
    params.anals.chans = [1:30];
    params.anals.winsecsGAMMA = 4; 
    [params.anals.nozFilter_bGAMMA,params.anals.nozFilter_aGAMMA] = nozFilter(1,1.5); 
    params.preproc.refType = 'tps'; %tps-TP9 + TP10; avg average ref
elseif strcmp(preproc,'p5')
    %same preproc as 'p2'
    params.anals.winsecs = 120; %used for theta and delta
    [b,a] = nozFilter(1,3); %Filter for local gain calc, Nozaraden-like
    params.anals.chans = [1:32];
    params.anals.winsecsGAMMA = 4; 
    [params.anals.nozFilter_bGAMMA,params.anals.nozFilter_aGAMMA] = nozFilter(1,1.5); 
    params.preproc.refType = 'avg'; %tps-TP9 + TP10; avg average ref
end


%%
params.chanNames = ["Fp1","Fz","F3","F7","FT9","FC5","FC1","C3","T7","TP9","CP5","CP1","Pz","P3","P7","O1","Oz","O2","P4","P8","TP10","CP6","CP2","Cz","C4","T8","FT10","FC6","FC2","F4","F8","Fp2"]; %what jakobs is
params.trigOffset = {'S102' 'S104' 'S002' 'S004'}; %doesn't have to be in chron order same for everyone
params.chanRange = 1:32;
%params.refChannel = 'A2'; not sure about this
params.gammaSongs = {'GammaBoardwalk';'GammaGraceland'};
params.regSongs = {'RegBoardwalk';'RegGraceland'};

params.anals.chans2use = chans2use; %channels to us for PSD
%params.anals.winsecs = 120; %window size for PSD 
%params.anals.chans = [1:30]; %
%[b,a] = nozFilter(1,3); %Filter for local gain calc, Nozaraden-like 
params.anals.nozFilter_a = a;
params.anals.nozFilter_b = b;
%subject specific gamma frequencies setGamma = box fixed @ 39 hz after 10/25
stimBasedGammaBoardwalk = 40.968;
stimBasedGammaGraceland = 39.2150;
%{
frequencies from light timesereis peaks delta*20
stimBasedGammaBoardwalk = 40.2832;
stimBasedGammaGraceland = 39.0625;
stimBasedGammaBoardwalk = 36.926270; %41.015625;
stimBasedGammaGraceland = 37.04834;
%}

setGamma = 39.000;%38.610;%39.001465;
params.anals.freqsOfInterest = [2.0484, 4.0968, 19.500, setGamma;          % boardwalk Stimulus frequencies
                                 1.96075, 3.9215, 19.500, setGamma;];       % graceland Stimulus frequencies
 %{
%frequencies from light timesereis peaks
params.anals.freqsOfInterest = [2.01416, 4.089355, 19.500, setGamma;          % boardwalk Stimulus frequencies
                                 1.953125, 3.967285, 19.500, setGamma;];       % graceland Stimulus frequencies
 %}

params.anals.flims = [1 50];                         % Frequency limits for PSD plots
params.anals.clims = 1.2*[-1 1];                     % Gain limits for topos and local gain orig 1.5*

%%%%%%%% Remove raw EEG before experiment, between trials, and after experiment? %%%%%%
% 'stimMarkers' = pass through stimulus codes
% 'stimPnts' = pass through sample points
% 'false' = don't epoch
params.preproc.initialEpoch = 'stimMarkers';
%%%%%% Epoch parameters(May not be used) %%%%%%
params.preproc.prepostStim = -0.5; %was just preStim orig, but used in code to extend both pre and post EEG epoch by amount (should be neg)
%%%%%% Filter parameters %%%%%%
params.preproc.lpCutoff = 55;        % Cut-off freq (Hz)for low-pass filter, set to nan if dn't want
params.preproc.hpCutoff = 1;         % Cut-off freq (Hz) for high-pass filter, set to nan if dn't want
params.preproc.notchFreq = 60; 
params.preproc.downSampleTo = 500;
%params.preproc.refType = 'avg'; %tps TP9 + TP10; avg average ref

%% subject info parameters (subjects to run)
isub=0;

params.preproc.setting = preproc;

%{
isub=isub+1;
subjectInfo(isub).batch = '';
subjectInfo(isub).subjectID = '';
subjectInfo(isub).intervention = "";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{};{}}; %in chron order
subjectInfo(isub).preciseTriggers = ; %true or false
subjectInfo(isub).origRawFolders = {'';''};
subjectInfo(isub).origRawFiles = {'';''};
subjectInfo(isub).rawFiles2use = {sprintf('%sevntFXD',preproc);sprintf('%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {,; , }; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-10;15}; %used in plotOrganize.m
%}

%% pilot (RA) subs
% Aug 23 Glasses MacGyver Pilot Batch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
isub=isub+1;
subjectInfo(isub).batch = 'MacGyver';
subjectInfo(isub).subjectID = 'EWU';
subjectInfo(isub).intervention = "NA";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).origRawFolders = {'230816EWU'};
subjectInfo(isub).origRawFiles = {'230816EWU_glasses_prototype'};
subjectInfo(isub).rawFiles2use = {sprintf('230816EWU_glasses_prototype%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230816EWU_glasses_prototype%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma; setGamma }; 
subjectInfo(isub).powerYlims = {-15;15}; %used in plotOrganize.m
%}

% EVT1 Glasses Pilot Batch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
isub=isub+1;
subjectInfo(isub).batch = 'EVT1';
subjectInfo(isub).subjectID = 'PLOU';
subjectInfo(isub).intervention = "NA";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).origRawFolders = {'230816EWU'};
subjectInfo(isub).origRawFiles = {'230816EWU_glasses_prototype'};
subjectInfo(isub).rawFiles2use = {sprintf('230816EWU_glasses_prototype%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230816EWU_glasses_prototype%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma; setGamma }; 
subjectInfo(isub).powerYlims = {-15;15}; %used in plotOrganize.m
%}

% RA DRY Pilot Batch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
isub=isub+1;
subjectInfo(isub).batch = 'DRYRA2024';
subjectInfo(isub).subjectID = 'ANNI';
subjectInfo(isub).intervention = "NA";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre'};
subjectInfo(isub).allTrig = {{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).origRawFolders = {'240906ANNI1'};
subjectInfo(isub).origRawFiles = {'240906ANNI1'};
subjectInfo(isub).rawFiles2use = {sprintf('240906ANNI1%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240906ANNI1%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma; setGamma }; 
subjectInfo(isub).powerYlims = {-25;-5}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'DRYRA2024';
subjectInfo(isub).subjectID = 'QNCY';
subjectInfo(isub).intervention = "NA";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).origRawFolders = {'240913QNCY1'};
subjectInfo(isub).origRawFiles = {'240913QNCY1'};
subjectInfo(isub).rawFiles2use = {sprintf('240913QNCY1%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240913QNCY1%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma; setGamma }; 
subjectInfo(isub).powerYlims = {-25;-5}; %used in plotOrganize.m
%}

% RA Pilot Batch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{ 
isub=isub+1;
subjectInfo(isub).batch = 'RA2024';
subjectInfo(isub).subjectID = 'ANNI';
subjectInfo(isub).intervention = "NA";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre'};
subjectInfo(isub).allTrig = {{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).origRawFolders = {'240410ANNI1'};
subjectInfo(isub).origRawFiles = {'240410ANNI1'};
subjectInfo(isub).rawFiles2use = {sprintf('240410ANNI1%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240410ANNI1%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma; setGamma }; 
subjectInfo(isub).powerYlims = {-25;-5}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'RA2024';
subjectInfo(isub).subjectID = 'QNCY';
subjectInfo(isub).intervention = "NA";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).origRawFolders = {'240210QNCY3'};
subjectInfo(isub).origRawFiles = {'240410QNCY3'};
subjectInfo(isub).rawFiles2use = {sprintf('240410QNCY3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240410QNCY3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma; setGamma }; 
subjectInfo(isub).powerYlims = {-25;-5}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'RA2024';
subjectInfo(isub).subjectID = 'JNYU';
subjectInfo(isub).intervention = "NA";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).origRawFolders = {'240412JNYU1'};
subjectInfo(isub).origRawFiles = {'240412JNYU1'};
subjectInfo(isub).rawFiles2use = {sprintf('240412JNYU1%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240412JNYU1%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma; setGamma }; 
subjectInfo(isub).powerYlims = {-15;15}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'RA2024';
subjectInfo(isub).subjectID = 'BOWN';
subjectInfo(isub).intervention = "NA";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre'};
subjectInfo(isub).allTrig = {{'S003' 'S004' 'S001' 'S002' 'S101' 'S102' 'S103' 'S104'}}; %in chron order
subjectInfo(isub).origRawFolders = {'240410BOWN3'};
subjectInfo(isub).origRawFiles = {'240410BOWN1'};
subjectInfo(isub).rawFiles2use = {sprintf('240410BOWN1%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240410BOWN1%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma; setGamma }; 
subjectInfo(isub).powerYlims = {-15;15}; %used in plotOrganize.m
%}

%% REAL subs
% "Gamma MBI Newbridge Batch" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
isub=isub+1;
subjectInfo(isub).batch = 'Newbridge'; %Newbridge
subjectInfo(isub).subjectID = 'PWOL';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'};{'S001' 'S002' 'S003' 'S004' 'S101' 'S102' 'S103' 'S104'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'230324PWOL1';'230606PWOL3'};
subjectInfo(isub).origRawFiles = {'230324PWOL1';'230606PWOL3'};
subjectInfo(isub).rawFiles2use = {sprintf('230324PWOL1%sevntFXD',preproc);sprintf('230606PWOL3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230324PWOL1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230606PWOL3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-15;15}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'Newbridge'; %Newbridge
subjectInfo(isub).subjectID = 'FWOL';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S001' 'S002' 'S003' 'S004' 'S101' 'S102' 'S103' 'S104'};{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'230327FWOL1';'230606FWOL3'};
subjectInfo(isub).origRawFiles = {'230327FWOL1';'230606FWOL3'};
subjectInfo(isub).rawFiles2use = {sprintf('230327FWOL1%sevntFXD',preproc);sprintf('230606FWOL3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230327FWOL1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230606FWOL3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-15;15}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'Newbridge'; %Newbridge
subjectInfo(isub).subjectID = 'SSCH';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'};{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
ssubjectInfo(isub).preciseTriggers = false;
ubjectInfo(isub).origRawFolders = {'230331SSCH1';'230620SSCH3'};
subjectInfo(isub).origRawFiles = {'230331SSCH1';'230620SSCH3'}; %VERIFY SESS 1 !!!!!!
subjectInfo(isub).rawFiles2use = {sprintf('230331SSCH1%sevntFXD',preproc);sprintf('230620SSCH3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230331SSCH1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230620SSCH3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-20;10}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'Newbridge'; %Newbridge
subjectInfo(isub).subjectID = 'BGOL';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S001' 'S002' 'S003' 'S004' 'S103' 'S104' 'S101' 'S102'};{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'230407BGOL1';'230623BGOL3'};
subjectInfo(isub).origRawFiles = {'230407BGOL1';'230623BGOL3'};
subjectInfo(isub).rawFiles2use = {sprintf('230407BGOL1%sevntFXD',preproc);sprintf('230623BGOL3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230407BGOL1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230623BGOL3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-10;15}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'Newbridge'; %Newbridge
subjectInfo(isub).subjectID = 'BROS';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'};{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'230411BROS1';'230718BROS3'};
subjectInfo(isub).origRawFiles = {'230411BROS1';'230718BROS3'};
subjectInfo(isub).rawFiles2use = {sprintf('230411BROS1%sevntFXD',preproc);sprintf('230718BROS3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230411BROS1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230718BROS3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-15;15}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'Newbridge'; %Newbridge
subjectInfo(isub).subjectID = 'SSHE';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S001' 'S002' 'S003' 'S004' 'S101' 'S102' 'S103' 'S104'};{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'230411SSHE1';'230718SSHE3'};
subjectInfo(isub).origRawFiles = {'230411SSHE1_try2';'230718SSHE3'}; %
subjectInfo(isub).rawFiles2use = {sprintf('230411SSHE1_try2%sevntFXD',preproc);sprintf('230718SSHE3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230411SSHE1_try2%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230718SSHE3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-15;20}; %used in plotOrganize.m
%}

%% start subjects who have been preprocessed
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "Gamma MBI First Batch" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2022
subjectInfo(isub).subjectID = 'DMUN';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'};{'S001' 'S002' 'S003' 'S004' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'220811DMUN1';'220811DMUN3'};
subjectInfo(isub).origRawFiles = {'220811DMUN1';'220811DMUN3'};
subjectInfo(isub).rawFiles2use = {sprintf('220811DMUN1%sevntFXD',preproc);sprintf('220811DMUN3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('220811DMUN1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('220811DMUN3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2022
subjectInfo(isub).subjectID = 'JDEN';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S003' 'S004' 'S001' 'S002'};{'S001' 'S002' 'S003' 'S004' 'S101' 'S102' 'S103' 'S104'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'220919JDEN1';'220919JDEN3'};
subjectInfo(isub).origRawFiles = {'220919JDEN1';'220919JDEN3'};
subjectInfo(isub).rawFiles2use = {sprintf('220919JDEN1%sevntFXD',preproc);sprintf('220919JDEN3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('220919JDEN1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('220919JDEN3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2022
subjectInfo(isub).subjectID = 'EKAS';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S001' 'S002' 'S003' 'S004'};{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'221103EKAS1';'221103EKAS3'};
subjectInfo(isub).origRawFiles = {'221103EKAS1';'221103EKAS3'};
subjectInfo(isub).rawFiles2use = {sprintf('221103EKAS1%sevntFXD',preproc);sprintf('221103EKAS3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('221103EKAS1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('221103EKAS3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2022
subjectInfo(isub).subjectID = 'JMCG';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S003','S004' 'S001' 'S002' 'S101' 'S102' 'S103' 'S104'};{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'221024JMCG1';'221024JMCG3'};
subjectInfo(isub).origRawFiles = {'221024JMCG1';'221024JMCG3'};
subjectInfo(isub).rawFiles2use = {sprintf('221024JMCG1%sevntFXD',preproc);sprintf('221024JMCG3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('221024JMCG1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('221024JMCG3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-30;-5}; %used in plotOrganize.m

%{
isub=isub+1; %has trouble running ICA (crazy slow) didin't run
subjectInfo(isub).batch = 'GLMBI'; %ISEC2022
subjectInfo(isub).subjectID = 'GHER';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S001' 'S002' 'S003','S004' 'S103' 'S104' 'S101','S102'};{'S101' 'S102' 'S103' 'S104' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'220914GHER1';'220914GHER3'};
subjectInfo(isub).origRawFiles = {'220914GHER1';'GHER3'};
subjectInfo(isub).rawFiles2use = {sprintf('220914GHER1%sevntFXD',preproc);sprintf('GHER3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('220914GHER1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('GHER3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,stimBasedGammaBoardwalk; stimBasedGammaGraceland,stimBasedGammaGraceland}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-15;15}; %used in plotOrganize.m
%}

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "Gamma MBI Most Recent Batch" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'JPEN';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'};{'S001' 'S002' 'S003' 'S004' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'231005JPEN1';'231005JPEN3'};
subjectInfo(isub).origRawFiles = {'231005JPEN1';'231005JPEN3'};
subjectInfo(isub).rawFiles2use = {sprintf('231005JPEN1%sevntFXD',preproc);sprintf('231005JPEN3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('231005JPEN1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('231005JPEN3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,setGamma; stimBasedGammaGraceland,setGamma }; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'PMUL';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S003' 'S004' 'S001' 'S002'};{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'230925PMUL1';'230925PMUL3'};
subjectInfo(isub).origRawFiles = {'230925PMUL1';'230925PMUL3'};
subjectInfo(isub).rawFiles2use = {sprintf('230925PMUL1%sevntFXD',preproc);sprintf('230925PMUL3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230925PMUL1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230925PMUL3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,setGamma; stimBasedGammaGraceland,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'LJOS';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'};{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'230918LJOS1';'230918LJOS3'};
subjectInfo(isub).origRawFiles = {'230918LJOS1';'230918LJOS3'};
subjectInfo(isub).rawFiles2use = {sprintf('230918LJOS1%sevntFXD',preproc);sprintf('230918LJOS3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230918LJOS1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230918LJOS3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,setGamma; stimBasedGammaGraceland,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'JPRI';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'};{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'230913JPRI1';'230913JPRI3'};
subjectInfo(isub).origRawFiles = {'230913JPRI1';'231211JPRI3'};
subjectInfo(isub).rawFiles2use = {sprintf('230913JPRI1%sevntFXD',preproc);sprintf('231211JPRI3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230913JPRI1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('231211JPRI3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,setGamma; stimBasedGammaGraceland,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-20;5}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'RFIS';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'};{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'230911RFIS1';'230911RFIS3'};
subjectInfo(isub).origRawFiles = {'230911RFIS1';'230911RFIS3'};
subjectInfo(isub).rawFiles2use = {sprintf('230911RFIS1%sevntFXD',preproc);sprintf('230911RFIS3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230911RFIS1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230911RFIS3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,setGamma; stimBasedGammaGraceland,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-30;-5}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'BYET';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'};{'S001' 'S002' 'S003' 'S004' 'S101' 'S102' 'S103' 'S104'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'230906BYET1';'230831BYET3'};
subjectInfo(isub).origRawFiles = {'230906BYET1';'231120BYET3'};
subjectInfo(isub).rawFiles2use = {sprintf('230906BYET1%sevntFXD',preproc);sprintf('231120BYET3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230906BYET1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('231120BYET3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,setGamma; stimBasedGammaGraceland,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'ACON';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S003' 'S004' 'S001' 'S002' 'S101' 'S102' 'S103' 'S104'};{'S101' 'S102' 'S103' 'S104' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'230907ACON1';'230906ACON3'};
subjectInfo(isub).origRawFiles = {'230907ACON1';'230906ACON3'};
subjectInfo(isub).rawFiles2use = {sprintf('230907ACON1%sevntFXD',preproc);sprintf('230906ACON3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230907ACON1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230906ACON3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,setGamma; stimBasedGammaGraceland,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-35;-10}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'SRE';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S001' 'S002' 'S003' 'S004' 'S103' 'S104' 'S101' 'S102'};{'S101' 'S102' 'S103' 'S104' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'230922SRE1';'230922SRE3'};
subjectInfo(isub).origRawFiles = {'230922SRE1';'230922SRE3'};
subjectInfo(isub).rawFiles2use = {sprintf('230922SRE1%sevntFXD',preproc);sprintf('230922SRE3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('230922SRE1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('230922SRE3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {stimBasedGammaBoardwalk,setGamma; stimBasedGammaGraceland,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-30;-5}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'JCHA';
subjectInfo(isub).intervention = "MBI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'};{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'231013JCHA1';'231013JCHA3'};
subjectInfo(isub).origRawFiles = {'231013JCHA1';'231013JCHA3'};
subjectInfo(isub).rawFiles2use = {sprintf('231013JCHA1%sevntFXD',preproc);sprintf('231013JCHA3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('231013JCHA1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('231013JCHA3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-35;-10}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'STER';
subjectInfo(isub).intervention = "MBI"; 
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S001' 'S002' 'S003' 'S004' 'S103' 'S104' 'S101' 'S102'};{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'240926STER1';'240926STER3'};
subjectInfo(isub).origRawFiles = {'260926STER1';'240926STER3'};
subjectInfo(isub).rawFiles2use = {sprintf('260926STER1%sevntFXD',preproc);sprintf('240926STER3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('260926STER1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240926STER3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m


isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'BTOR';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'};{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'231025BTOR1';'231025BTOR3'};
subjectInfo(isub).origRawFiles = {'231025BTOR1';'231025BTOR3'};
subjectInfo(isub).rawFiles2use = {sprintf('231025BTOR1%sevntFXD',preproc);sprintf('231025BTOR3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('231025BTOR1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('231025BTOR3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/25
subjectInfo(isub).powerYlims = {-20;-5}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2023
subjectInfo(isub).subjectID = 'SHEF';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S003' 'S004' 'S001' 'S002' 'S103' 'S104' 'S101' 'S102'};{'S003' 'S004' 'S001' 'S002' 'S101' 'S102' 'S103' 'S104'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'231115SHEF1';'231115SHEF3'};
subjectInfo(isub).origRawFiles = {'231115SHEF1_2';'bp-sequencer_001'};
subjectInfo(isub).rawFiles2use = {sprintf('231115SHEF1%s_2evntFXD',preproc);sprintf('bp-sequencer_001%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('231115SHEF1%s_2evntFXD_ResampFiltRejInterpRef',preproc);sprintf('bp-sequencer_001%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'LPYK';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'};{'S103' 'S104' 'S101' 'S102' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).preciseTriggers = false;
subjectInfo(isub).origRawFolders = {'240118LPYK1';'240118LPYK3'};
subjectInfo(isub).origRawFiles = {'240118LPYK1';'240118LPYK3'};
subjectInfo(isub).rawFiles2use = {sprintf('240118LPYK1%sevntFXD',preproc);sprintf('240118LPYK3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240118LPYK1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240118LPYK3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'DPAS';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'};{'S001' 'S002' 'S003' 'S004' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'240306DPAS1';'240306DPAS3'};
subjectInfo(isub).origRawFiles = {'24036DPAS1';'240306DPAS3'};
subjectInfo(isub).rawFiles2use = {sprintf('24036DPAS1%sevntFXD',preproc);sprintf('240306DPAS3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('24036DPAS1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240306DPAS3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-35;-10}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'PMAN';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S001' 'S002' 'S003' 'S004' 'S103' 'S104' 'S101' 'S102'};{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'240328PMAN1';'240529PMAN3'};
subjectInfo(isub).origRawFiles = {'240328PMAN1';'240328PMAN3'};
subjectInfo(isub).rawFiles2use = {sprintf('240328PMAN1%sevntFXD',preproc);sprintf('240328PMAN3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240328PMAN1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240328PMAN3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'BFIT';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S001' 'S002' 'S003' 'S004'};{'S003' 'S004' 'S001' 'S002' 'S101' 'S102' 'S103' 'S104' }}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'240410BFIT1';'240410BFIT3'};
subjectInfo(isub).origRawFiles = {'240410BFIT1';'240410BFIT3'};
subjectInfo(isub).rawFiles2use = {sprintf('240410BFIT1%sevntFXD',preproc);sprintf('240410BFIT3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240410BFIT1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240410BFIT3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-20;-5}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'LPAU';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S003' 'S004' 'S001' 'S002'};{'S103' 'S104' 'S101' 'S102' 'S001' 'S002' 'S003' 'S004'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'240516LPAU1';'240516LPAU3'};
subjectInfo(isub).origRawFiles = {'240516LPAU1';'240516LPAU3'};
subjectInfo(isub).rawFiles2use = {sprintf('240516LPAU1%sevntFXD',preproc);sprintf('240516LPAU3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240516LPAU1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240516LPAU3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-35;-10}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'MHUA';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S001' 'S002' 'S003' 'S004'};{'S003' 'S004' 'S001' 'S002','S101' 'S102' 'S103' 'S104'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'240828MHUA1';'240828MHUA3'};
subjectInfo(isub).origRawFiles = {'240828MHUA1';'240828MHUA3'};
subjectInfo(isub).rawFiles2use = {sprintf('240828MHUA1%sevntFXD',preproc);sprintf('240828MHUA3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240828MHUA1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240828MHUA3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'DWEL';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S001' 'S002' 'S003' 'S004' 'S101' 'S102' 'S103' 'S104'};{'S001' 'S002' 'S003' 'S004' 'S103' 'S104' 'S101' 'S102'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'240904DWEL1';'240904DWEL3'};
subjectInfo(isub).origRawFiles = {'240904DWEL1';'240904DWEL3'};
subjectInfo(isub).rawFiles2use = {sprintf('240904DWEL1%sevntFXD',preproc);sprintf('240904DWEL3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240904DWEL1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240904DWEL3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'BROB';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S001' 'S002' 'S003' 'S004'};{'S001' 'S002' 'S003' 'S004' 'S101' 'S102' 'S103' 'S104'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'240925BROB1';'240925BROB3'};
subjectInfo(isub).origRawFiles = {'240925BROB1';'240925BROB3'};
subjectInfo(isub).rawFiles2use = {sprintf('240925BROB1%sevntFXD',preproc);sprintf('240925BROB3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('240925BROB1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('240925BROB3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'SHUF';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S103' 'S104' 'S101' 'S102' 'S001' 'S002' 'S003' 'S004'};{'S001' 'S002' 'S003' 'S004' 'S101' 'S102' 'S103' 'S104'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'241115SHUF1';'241115SHUF3'};
subjectInfo(isub).origRawFiles = {'241115SHUF1';'241115SHUF3'};
subjectInfo(isub).rawFiles2use = {sprintf('241115SHUF1%sevntFXD',preproc);sprintf('241115SHUF3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('241115SHUF1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('241115SHUF3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

isub=isub+1;
subjectInfo(isub).batch = 'GLMBI'; %ISEC2024
subjectInfo(isub).subjectID = 'JPAT';
subjectInfo(isub).intervention = "HEI";
subjectInfo(isub).swapBoxes = 'no'; %'yes' if EEG boxes were plugged in wrong, 'no' otherwise
subjectInfo(isub).sessNames = {'Pre';'Post'};
subjectInfo(isub).allTrig = {{'S101' 'S102' 'S103' 'S104' 'S003' 'S004' 'S001' 'S002'};{'S103' 'S104' 'S101' 'S102' 'S003' 'S004' 'S001' 'S002'}}; %in chron order
subjectInfo(isub).preciseTriggers = true;
subjectInfo(isub).origRawFolders = {'241120JPAT1';'241120JPAT3'};
subjectInfo(isub).origRawFiles = {'241120JPAT1';'241120JPAT3'};
subjectInfo(isub).rawFiles2use = {sprintf('241120JPAT1%sevntFXD',preproc);sprintf('241120JPAT3%sevntFXD',preproc)}; %post correcting events
subjectInfo(isub).files4ica = {sprintf('241120JPAT1%sevntFXD_ResampFiltRejInterpRef',preproc);sprintf('241120JPAT3%sevntFXD_ResampFiltRejInterpRef',preproc)}; %post preproc_timeseries.m
subjectInfo(isub).gammFreqs_B_G = {setGamma,setGamma; setGamma,setGamma}; %pre and post %diff box fixed @ 39 after 10/26
subjectInfo(isub).powerYlims = {-25;0}; %used in plotOrganize.m

%% end subjects who have been preprocessed


%%
function [trigCode] = makeStimulusCodes(makeTrig, omitTrig)
% Parker Tichko, 2020
% makeTrig is a scalar or vector of integers to make stimulus codes for
% omitTrig is a scalar of vector of integer to omit
% Example: makeStimulusCodes(1, 0) returns "S  1"
% Example: makeStimulusCodes(1:3, 0) returns "S  1", "S  2", "S  3"
% Example: makeStimulusCodes(1:3, 2) returns "S  1", "S  3"

trigCode = [];

    for i = 1:length(makeTrig) 
        if makeTrig(i) ~= omitTrig
            if strlength(int2str(makeTrig(i))) == 1
                trigCode = [trigCode {['S  ', int2str(makeTrig(i))]}];
            elseif strlength(int2str(makeTrig(i))) == 2
                  trigCode = [trigCode {['S ', int2str(makeTrig(i))]}];
            elseif strlength(int2str(makeTrig(i))) == 3
                  trigCode = [trigCode {['S', int2str(makeTrig(i))]}];
            end
        end
    end
    
end

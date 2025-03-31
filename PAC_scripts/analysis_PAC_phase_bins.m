function [outData] = analysis_PAC_phase_bins(ringCond,interventionType,batch,chans2use,subjectInfo,params)
 %{
    AA - Nov 2024
 %}

if length(chans2use)~=0
    chnstr = sprintf('%dchns',length(chans2use));
else
    chnstr = 'allchns'; %using all
end
%init table to hold some analysis info.
varTypes = ["string","string","string","double"]; %repelem(["double"],32)];
varNames = ["subjectID","fileName","Session","TotalTime"];%params.chanNames];

outTable = table('Size',[length({subjectInfo.subjectID})*2,length(varTypes)],'VariableTypes',varTypes,'VariableNames',varNames);
dumbcounter = 1;

% PAC PARAMS:
amplitudeFreqRange = [10 50]; numAmpFreqs = 41;

num_phase_bins = 18;
phase_bins = linspace(-pi,pi, num_phase_bins+1);

surrogate_size = 200;

method_complex = 'wavelet';
num_cycles = 5;
surrogate = true;

highFreqs = linspace(amplitudeFreqRange(1), amplitudeFreqRange(2), numAmpFreqs);

%assuming 30 channels
% alldata = zeros(length(subjectInfo(1).sessNames), length(params.gammaSongs), length({subjectInfo.subjectID}),30, numPhaseFreqs, numAmpFreqs);

numSubjects = length({subjectInfo.subjectID});

for isess=1:length(subjectInfo(1).sessNames) %length(subjectInfo(isub).origRawFiles)
    gmTopoD = [];
    gmTopoT = [];
    gmTopoG = [];

    if strcmp(ringCond,'GAM')
        songs = params.gammaSongs;
    end 
    if strcmp(ringCond,'REG')
        songs = params.regSongs;
    end 
    for isong = 1:length(songs) %songs)
        mLP  = [];              %  accumlating all elctrode average Power
        mG   = [];              %  accumlating all elctrode average Gain
        mTopoD = [];            %  mean topos over song Gain
        mTopoT = [];
        mTopoG = [];
        submeanLinPs = []; % list of subject song linear psds for use in fooof
        fprintf('working on %s...\n',songs{isong})
        song = songs{isong};
        if contains(song, "Boardwalk")
            phaseFreqs = [params.anals.freqsOfInterest(1, 2) 2*params.anals.freqsOfInterest(1, 2)];
        elseif contains(song, "Graceland")
            phaseFreqs = [params.anals.freqsOfInterest(2, 2) 2*params.anals.freqsOfInterest(2, 2)];
        else
            error("song not recognized");
        end
        
        for isub = 1:numSubjects %subjs
            subj = subjectInfo(isub).subjectID;
            
            song = songs{isong};

            fprintf('working on session: %s, subject: %s\n',subjectInfo(isub).sessNames{isess},subjectInfo(isub).subjectID)

            tmpnf_eeg = dir(fullfile(params.path2data,subjectInfo(isub).sessNames{isess},'preprocessed',subjectInfo(isub).origRawFolders{isess},[subjectInfo(isub).files4ica{isess} 'Ica' song '.set']));

            if ~isempty({tmpnf_eeg.name}) && length({tmpnf_eeg.name})==1 && subjectInfo(isub).intervention==interventionType

                %% setup
                %make sure preprocessed dir exists for subject
                savepath = fullfile(params.path2data,subjectInfo(isub).sessNames{isess},'analyzed',subjectInfo(isub).origRawFolders{isess}); %[preprocEEGPath participantID];
                mkdir(savepath);

                EEG = pop_loadset('filename', tmpnf_eeg.name,'filepath',tmpnf_eeg.folder);
                
                if length(chans2use)~=0
                    chnidxs = [];
                    for ichn=1:length(chans2use)
                        chnidxs(end+1) = find(strcmpi({EEG.chanlocs.labels},chans2use(ichn))==1);
                    end
                    EEG = pop_select(EEG,'channel',chnidxs);
                else
                    chnidxs = params.anals.chans;

                end

                labels = {EEG.chanlocs.labels};
                win = EEG.srate*params.anals.winsecs;
                ovr = win/2;

                % X  = EEG.data;
                % disp("made it here")

               
                
                disp("computing coupling now: ")
                
                
                eegdata = EEG.data';
                eeg_srate = EEG.srate;  

                final_amplitudes_by_phases = zeros(30,length(phaseFreqs), num_phase_bins, length(highFreqs));
                raw_amplitudes_by_phases = zeros(30,length(phaseFreqs), num_phase_bins, length(highFreqs));

                tic
                parfor low_idx = 1:length(phaseFreqs)
                    disp(strcat("   low freq: ", num2str(phaseFreqs(low_idx))));
 

                    [low_freq_phases, ~] = eegConvolution(eegdata, num_cycles, phaseFreqs(low_idx), eeg_srate);
                    
                    partition_lims = [round(1/6 * size(low_freq_phases, 1)) round(5/6 * size(low_freq_phases, 1))]; % make sure there is a sufficiently large lag

                    for high_idx = 1:numAmpFreqs
                        disp(highFreqs(high_idx));
                        [~, all_amplitudes] = eegConvolution(eegdata, num_cycles, highFreqs(high_idx), eeg_srate);
                        surrogate_amplitudes_by_phases = zeros(30, surrogate_size, num_phase_bins);
                        amplitudes_by_phases = zeros(30, num_phase_bins);
                        if surrogate
                            surrogate_amplitudes = zeros(30, surrogate_size, size(all_amplitudes, 1));
                            for i = 1:surrogate_size
                                surrogate_amplitudes(:,i, :) = circshift(all_amplitudes, randi(partition_lims), 1)';
                            end
                        end
                        for bin = 1:num_phase_bins
                            for electrode = 1:30
                                amplitudes_by_phases(electrode, bin) = mean(all_amplitudes(low_freq_phases(:,electrode) > phase_bins(bin) & ...
                                    low_freq_phases(:,electrode) < phase_bins(bin+1), electrode), 1);
                                if surrogate
                                    surrogate_amplitudes_by_phases(electrode,:, bin) = mean(surrogate_amplitudes(electrode,:,low_freq_phases(:,electrode) > phase_bins(bin) & ...
                                        low_freq_phases(:,electrode) < phase_bins(bin+1)), 3);
                                end
                            end
                        end
                        raw_amplitudes_by_phases(:,low_idx,:,high_idx) = amplitudes_by_phases; 
                        if surrogate
                            mean_surrogate = squeeze(mean(surrogate_amplitudes_by_phases, 2));
                            std_surrogate = squeeze(std(surrogate_amplitudes_by_phases, [],2));
                            final_amplitudes_by_phases(:,low_idx,:,high_idx) = (amplitudes_by_phases - mean_surrogate)./std_surrogate;                        
                        end
                    end
   
                end
                toc

                save_couplings(savepath, strcat('PAC_phase_bins_',subjectInfo(isub).sessNames{isess},'_', songs{isong},'_',subjectInfo(isub).subjectID), raw_amplitudes_by_phases, final_amplitudes_by_phases, phase_bins, amplitudeFreqRange, numAmpFreqs, phaseFreqs)

            end
        end
    end
end

%outData = alldata;
outData = [];
end

function save_couplings(savepath, filename, raw_amplitudes_by_phases, final_amplitudes_by_phases, phase_bins, amplitudeFreqRange, numAmpFreqs, phaseFreqs)
    save(fullfile(savepath, filename), "raw_amplitudes_by_phases", "final_amplitudes_by_phases", "phase_bins", "amplitudeFreqRange", "numAmpFreqs", "phaseFreqs")
end


% expects each column to be a channel timeseries 
function [phase_timeSeries, amplitude_timeSeries] = eegConvolution(eegdata, cycles, frequency, srate)
    time  = -3:1/srate:3; % best practice is to have time=0 at the center of the wavelet
    % create complex sine wave
    sine_wave = exp( 1i*2*pi*frequency.*time );
    % create Gaussian window
    s = cycles / (2*pi*frequency); % this is the standard deviation of the gaussian
    gaus_win  = exp( (-time.^2) ./ (2*s^2) );
    % now create Morlet wavelet
    cmw = gaus_win.*sine_wave;

    nData = length(eegdata);
    nKern = length(cmw);
    nConv = nData + nKern - 1;
    %FFTs:
    % note that the "N" parameter is the length of convolution, NOT the length
    % of the original signals! Super-important!
    % FFT of wavelet, and amplitude-normalize in the frequency domain
    cmwX = fft(cmw,nConv);
    cmwX = cmwX ./ max(cmwX);
    % FFT of data
    dataX = fft(eegdata,nConv);
    conv_res = dataX.*cmwX.'; %% --> we have just performed a convolution in the frequency domain (potentially the coolest thing ever) 
    % now back to the time domain
    % cut 1/2 of the length of the wavelet from the beginning and from the end
    half_wav = floor( length(cmw)/2 )+1;
    % take inverse Fourier transform
    conv_res_timedomain = ifft(conv_res);
    conv_res_timedomain = conv_res_timedomain(half_wav-1:end-half_wav,:);
    % conv_res_timedomain = conv_res_timedomain(100:end-100,:); %remove edge artifacts from convolution by trimming 100 samples from both ends
    phase_timeSeries = angle(conv_res_timedomain);
    amplitude_timeSeries = abs(conv_res_timedomain);
    % filtered_timeSeries = real(conv_res_timedomain);
    % power_timeSeries = abs(conv_res_timedomain).^2;
end

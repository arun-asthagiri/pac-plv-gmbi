function [outData] = analysis_window_PAC_fois(ringCond,interventionType,batch,chans2use,subjectInfo,params)
 %{
    AA - Jun 2024
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

% SLIDING PAC PARAMS:
cycles = @(freq) freq/5+3;
% take middle n seconds from song
excerpt_length = 120;
% sliding window
window_length_secs = 10;
hop_size_secs = 2;

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

        for isub = 1:numSubjects %subjs
            subj = subjectInfo(isub).subjectID;
            
            song = songs{isong};
            if contains(song, "Boardwalk")
                phaseFreqs = [params.anals.freqsOfInterest(1, 1) params.anals.freqsOfInterest(1, 2) 2*params.anals.freqsOfInterest(1, 2)];
            elseif contains(song, "Graceland")
                phaseFreqs = [params.anals.freqsOfInterest(2, 1) params.anals.freqsOfInterest(2, 2) 2*params.anals.freqsOfInterest(2, 2)];
            else
                error("song not recognized");
            end

            highFreq = subjectInfo(isub).gammFreqs_B_G{isong, isess};
        
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

                disp("computing coupling now: ")
                eegdata = EEG.data';
                eeg_srate = EEG.srate;  
                
                % identify indeces for middle n seconds
                idx_middle = floor(length(eegdata)/2);
                idx_start = idx_middle - (excerpt_length/2*eeg_srate) + 1;
                idx_end = idx_middle + (excerpt_length/2*eeg_srate);
                eegdata = eegdata(idx_start:idx_end,:);
                
                window_length = window_length_secs*eeg_srate;
                hop_size = hop_size_secs*eeg_srate;
                window_idxs = 1:hop_size:length(eegdata)-window_length-1;
                
                couplings = zeros(length(window_idxs), 30,length(phaseFreqs));
                raw_couplings = zeros(length(window_idxs), 30,length(phaseFreqs));
                tic
                for window_idx = 1:length(window_idxs) % iterate over windows
                    series = eegdata(window_idxs(window_idx):window_idxs(window_idx)+window_length, :);
                    [~, amplitude_high_freq] = eegConvolution(series, cycles(highFreq), highFreq, eeg_srate); % single high frequency
                    for low_idx = 1:length(phaseFreqs) % iterate over phase frequencies
                        disp(strcat("   low freq: ", num2str(phaseFreqs(low_idx))));
                        [low_freq_phases, ~] = eegConvolution(series, cycles(phaseFreqs(low_idx)), phaseFreqs(low_idx), eeg_srate);
                        aggregate_complex = amplitude_high_freq.* exp(1i*low_freq_phases);
                        coupling_vector = mean(aggregate_complex, 1);
                        [surrogate_mean, surrogate_std] = surrogatePAC(low_freq_phases, amplitude_high_freq, 200);
                        couplings(window_idx,:,low_idx) = (abs(coupling_vector)-surrogate_mean)./surrogate_std;
                        raw_couplings(window_idx,:,low_idx) = abs(coupling_vector);
                    end
                end
                toc

                save_couplings(savepath, strcat('PAC_sliding_fois_',subjectInfo(isub).sessNames{isess},'_', songs{isong},'_',subjectInfo(isub).subjectID), couplings, raw_couplings, phaseFreqs, highFreq,excerpt_length,window_length,hop_size, interventionType)
            

            end
        end
    end
end

%outData = alldata;
outData = [];
end

function save_couplings(savepath, filename, couplings, raw_couplings, phaseFreqs, highFreq, excerpt_length,window_length,hop_size, interventionType)
    save(fullfile(savepath, filename), "couplings", "raw_couplings", "phaseFreqs", "highFreq", "excerpt_length", "window_length", "hop_size", "interventionType")
end


function [mean_surrogate, std_surrogate] = surrogatePAC(phase_series, amplitude_series, num_permutations)
    partition_lims = [round(1/6 * length(phase_series)) round(5/6 * length(phase_series))]; % make sure there is a sufficiently large lag
    couplings = zeros(num_permutations, size(phase_series,2));

    for i = 1:num_permutations
        aggregate_complex = amplitude_series.* exp(1i*circshift(phase_series, randi(partition_lims), 1));
        couplings(i,:) = abs(mean(aggregate_complex, 1));
    end
    [mean_surrogate, std_surrogate] = deal(mean(couplings,1), std(couplings,1));
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

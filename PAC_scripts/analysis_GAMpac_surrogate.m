function [outData] = analysis_GAMpac_surrogate(ringCond,interventionType,batch,chans2use,subjectInfo,params)
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

% PAC PARAMS:
phaseFreqRange = [1 29]; numPhaseFreqs = 29;
amplitudeFreqRange = [30 54]; numAmpFreqs = 13;

method_complex = 'wavelet';
cycles = @(freq) freq/5+3;

surrogate = true;
lowFreqs = parallel.pool.Constant(linspace(phaseFreqRange(1), phaseFreqRange(2), numPhaseFreqs));
highFreqs = parallel.pool.Constant(linspace(amplitudeFreqRange(1), amplitudeFreqRange(2), numAmpFreqs));
% lowFreqs = linspace(phaseFreqRange(1), phaseFreqRange(2), numPhaseFreqs);
% highFreqs = linspace(amplitudeFreqRange(1), amplitudeFreqRange(2), numAmpFreqs);


% low_bandwidth = 1;
% high_bandwidth = 2;

if surrogate
    disp("running Surrogate with method: " + string(method_complex))
else
    disp("running without surrogate");
    warning("delete all parfor loops")
end

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

               
                if EEG.nbchan ~= 30
                    error("unexpected num chans")
                end

                disp("computing coupling now: ")
                
                
                eegdata = EEG.data';
                eegdata = parallel.pool.Constant(eegdata);
                eeg_srate = EEG.srate;  

                couplings = zeros(30,numPhaseFreqs, numAmpFreqs);
                raw_couplings = zeros(30,numPhaseFreqs, numAmpFreqs);
                peak_angles = zeros(30,numPhaseFreqs, numAmpFreqs);

                tic
                parfor low_idx = 1:numPhaseFreqs
                    disp(strcat("   low freq: ", num2str(lowFreqs.Value(low_idx))));

                    %if strcmp(method_complex, 'wavelet')
                    [low_freq_phases, ~] = eegConvolution(eegdata.Value, cycles(lowFreqs.Value(low_idx)), lowFreqs.Value(low_idx), eeg_srate);
                    %else
                     %   low_freq = bandpass(eegdata, [lowFreqs.Value(low_idx) - low_bandwidth/2, lowFreqs.Value(low_idx) + low_bandwidth/2 ], eeg_srate);
                      %  low_freq_phases = angle(hilbert(low_freq));
                    
                    %end
                    
                    
                    for high_idx = 1:numAmpFreqs
                        %if strcmp(method_complex, 'wavelet')
                        [~, amplitude_high_freq] = eegConvolution(eegdata.Value, cycles(highFreqs.Value(high_idx)), highFreqs.Value(high_idx), eeg_srate);
                        %else
                            %high_freq = bandpass(eegdata, [highFreqs.Value(high_idx)-high_bandwidth/2, highFreqs.Value(high_idx)+high_bandwidth/2], eeg_srate);
                            %amplitude_high_freq = abs(hilbert(high_freq));
                        %end
                        aggregate_complex = amplitude_high_freq.* exp(1i*low_freq_phases);
                        coupling_vector = mean(aggregate_complex, 1);
                        if surrogate
                            [surrogate_mean, surrogate_std] = surrogatePAC(low_freq_phases, amplitude_high_freq, 200);
                            couplings(:,low_idx, high_idx) = (abs(coupling_vector)-surrogate_mean)./surrogate_std;
                            raw_couplings(:,low_idx, high_idx) = abs(coupling_vector);
                        else
                            couplings(:,low_idx, high_idx) = abs(coupling_vector);
                        end
                        
                        peak_angles(:,low_idx,high_idx) = angle(coupling_vector);
                    end
                end
                toc

                save_couplings(savepath, strcat('PAC_Surrogate_cyc_change_',subjectInfo(isub).sessNames{isess},'_', songs{isong},'_',subjectInfo(isub).subjectID), couplings, raw_couplings, peak_angles, phaseFreqRange, numPhaseFreqs, amplitudeFreqRange, numAmpFreqs, interventionType)
            

            end
        end
    end
end

%outData = alldata;
outData = [];
end

function save_couplings(savepath, filename, couplings, raw_couplings,peak_angles, phaseFreqRange, numPhaseFreqs, amplitudeFreqRange, numAmpFreqs, interventionType)
    save(fullfile(savepath, filename), "couplings", "raw_couplings","peak_angles", "phaseFreqRange", "numPhaseFreqs", "amplitudeFreqRange", "numAmpFreqs", "interventionType")
end

function save_fig(fig, savepath, filename)
    saveas(fig, fullfile(savepath, filename));
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

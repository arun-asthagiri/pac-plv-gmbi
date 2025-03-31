function [outData] = analysis_PLV(ringCond,interventionType,batch,chans2use,subjectInfo,params)
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

boardwalk_path = '/work/mindlab/Projects/GammaMBI/Songs/Boardwalk.wav';
graceland_path = '/work/mindlab/Projects/GammaMBI/Songs/Graceland.wav';

% PAC PARAMS:
freqrange = [.2 15.2]; numFreqs = 76;
num_cycles = 5;
freqs = linspace(freqrange(1), freqrange(2), numFreqs);


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

            if contains(song, "Boardwalk")
                audio = miraudio(boardwalk_path);
                theta_freq = params.anals.freqsOfInterest(1, 2);
            elseif contains(song, "Graceland")
                audio = miraudio(graceland_path);
                theta_freq = params.anals.freqsOfInterest(2, 2);
            else
                error("song not recognized");
            end            
            % cochlear filtering
            resample_to = 500;
            filt = mirfilterbank(audio, 'Gammatone'); env = mirenvelope(filt, 'Sampling', resample_to);
            sum_bands = mirsum(env);
            env = get(sum_bands,'Data');
            audio_env = env{1}{1};
            audio_env_sr = get(sum_bands,'Sampling'); audio_env_sr = audio_env_sr{1};

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


                disp("computing coupling now: ")
                
                
                eegdata = EEG.data';
                eeg_srate = EEG.srate;  
                try
                    if (length(eegdata) - length(audio_env)) > 10 % check if more than 10 samples off
                        error("size mismatch")
                    end
                catch
                    disp("skipping participant due to size mismatch")
                    continue
                end
                % trim data 
                eegdata = eegdata(1:min(length(audio_env), length(eegdata)), :);
                audio_env = audio_env(1:min(length(audio_env), length(eegdata)), :);

                plvs = zeros(length(freqs), 30);
                
                tic
                for freq_idx = 1:length(freqs)
                    freq = freqs(freq_idx);
                    disp("   freq: "+ num2str(freq))
                    eeg_phases = eegConvolution(eegdata, num_cycles , freq, eeg_srate);
                    audio_phases = eegConvolution(audio_env, num_cycles , freq, audio_env_sr);
                    plvs(freq_idx,:) = abs(mean(exp(1i*(eeg_phases-audio_phases)), 1));
                end
                % calc theta foi PLV
                eeg_phases = eegConvolution(eegdata, num_cycles , theta_freq, eeg_srate);
                audio_phases = eegConvolution(audio_env, num_cycles , theta_freq, audio_env_sr);
                theta_plv = abs(mean(exp(1i*(eeg_phases-audio_phases)), 1));

                toc
                save(fullfile(savepath, strcat('PLV_high_res_',subjectInfo(isub).sessNames{isess},'_', songs{isong},'_',subjectInfo(isub).subjectID)),"plvs", "freqs","theta_plv", "interventionType")
            
            end
        end
    end
end
%outData = alldata;
outData = [];
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

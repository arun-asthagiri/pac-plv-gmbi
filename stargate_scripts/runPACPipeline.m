 %% phase-amplitude coupling analyses for gamma MBI EEG data %%
 %{
    BMK - April 2024
    AA - June 2024

 %}

function runPACPipeline(slurmIN)
%%
params_gMBI; eeglab redraw; close;

subjectInfo = subjectInfo(slurmIN);%
if iscell(unique({subjectInfo.batch})) && length(unique({subjectInfo.batch}))~=1
    fprintf('QUITING, found more than 1 batch type\n ADJUST PARAMS FILE\n')
    return
elseif length(unique({subjectInfo.batch}))==1
    batch = cell2mat( unique({subjectInfo.batch}));
    fprintf('working on batch: %s\n',batch)
end

% parpool(feature('numcores')); %for parallel computing on discovery cluster! https://rc-docs.northeastern.edu/en/latest/software/systemwide/matlab.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MBI PARTICIPANTS
[outTable_GAM] = analysis_compile_data('GAM','MBI', batch, params.anals.chans2use,subjectInfo,params);
disp(strcat("completed MBI gamma conditions: ", subjectInfo.subjectID))

% save(fullfile(params.homedir,"allDataPAC_oneFreq_Gam"),"outTable_GAM");
[outTable_REG] = analysis_compile_data('REG','MBI',batch,params.anals.chans2use,subjectInfo,params);
disp(strcat("completed MBI regular conditions: ", subjectInfo.subjectID))
% save(fullfile(params.homedir,"allDataPAC_oneFreq_Reg"),"outTable_REG");

% HEI Participants
[outTable_GAM] = analysis_compile_data('GAM','HEI', batch, params.anals.chans2use,subjectInfo,params);
disp(strcat("completed HEI gamma conditions: ", subjectInfo.subjectID))
% save(fullfile(params.homedir,"allDataPAC_oneFreq_Gam"),"outTable_GAM");
[outTable_REG] = analysis_compile_data('REG','HEI',batch,params.anals.chans2use,subjectInfo,params);
disp(strcat("completed HEI regular conditions: ", subjectInfo.subjectID))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot setup
% numSubjects = 9;
% numPhaseFreqs = 12; numAmpFreqs = 18;
% phaseFreqRange = [1, 12];
% amplitudeFreqRange = [20, 54];
% 
% 
% % 2x2x2: [Pre Post], [Reg Gamma], [Boardwalk Graceland]
% % all_couplings = zeros(numSubjects, 2, 2,2, 30, numPhaseFreqs, numAmpFreqs);
% all_couplings(:,:,2,:,:,:,:) = reshape(analsG_outTable_GAM.couplings, numSubjects, 2, 2, 30, numPhaseFreqs, numAmpFreqs);
% all_couplings(:,:,1,:,:,:,:) = reshape(analsG_outTable_REG.couplings, numSubjects, 2, 2, 30, numPhaseFreqs, numAmpFreqs);
%% draw plot
% figure;
% 
% % % 
% lowFreqs = linspace(phaseFreqRange(1), phaseFreqRange(2), numPhaseFreqs);
% highFreqs = linspace(amplitudeFreqRange(1), amplitudeFreqRange(2), numAmpFreqs);
% 
% electrodes_toplot = 1:30;
% 
% subjects_toplot = 1:numSubjects;
% 
% % for i = 1:30
% titles = {"Reg Boardwalk", "Reg Graceland", "Gamma Boardwalk", "Gamma Graceland"};
% 
% conditions = [1, 1, 1;   1, 1, 2;   1, 2, 1;  1, 2, 2]; % pre
% conditions = [2, 1, 1;   2, 1, 2;   2, 2, 1;  2, 2, 2]; % post
% 
% figure;
% 
% t = tiledlayout('flow');
% t.TileSpacing = "tight";
% 
% fontsize(15, "points");
% 
% xlabel(t,"frequency phase [Hz]");
% ylabel(t,"frequency amplitude [Hz]");
% for i = 1:length(titles)    
%     nexttile;
%     couplings = squeeze(mean(all_couplings(subjects_toplot,conditions(i,1),conditions(i,2),conditions(i,3),electrodes_toplot,:,:),[1, 5]));
%     pcolor(couplings');
%     shading interp
%     % cb = colorbar;
%     clim([-.02 .02])
% 
% 
%     % ylabel(cb,'Coupling Strength (au)','Rotation',270);
%     % % % 
%     xticks(lowFreqs);
%     yticks(1:length(highFreqs));
%     xticklabels(lowFreqs)
%     yticklabels(highFreqs)    
%     ylim([7 18])
%     title(titles{i});
%     fontsize(15, "points");
% end
% 
% title(t,"PAC: Post");
% t.Title.FontWeight = 'bold';
% t.Title.FontSize = 20;
% 
% % end
% c = colorbar;
% c.Layout.Tile = "east";
% 
% set(gcf,'color','w');

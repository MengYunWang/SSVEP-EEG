%% Pre-processing EEG data
% The preprocessing consisted of two stages:
% Stage I: preICA and run ICA
%   preICA: Import data> Channel Location > Filter > resample the data>  
%           rereference the data > extract all epoches > reject bad trials
%   run ICA: runica > exclude artifacts with ADjust
%   Baseline correction

% Stage II: extract conditions

% If you do not konw which function you should use, 
% please run one of your data, and type EEG.history.
% created by M.-Y. Wang
% 05-09-2017

%% Stage1-preICA
clear all
clc
cd ('F:\face-random\Raw_data')

for subi = 101:121;
    EEG = pop_biosig([num2str(subi),'-random.bdf'], 'channels',1:64,'ref',47);
    EEG.setname = [num2str(subi),'_preICA'];
    EEG = pop_chanedit(EEG, 'lookup','D:\\Program Files\\MATLAB\\R2014a\\matlabtoolbox\\eeglab2019_0\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp');
    EEG =pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',45);
    EEG = pop_resample( EEG, 1000);
    EEG = pop_reref( EEG, []);
    EEG = pop_epoch( EEG, {'1' '2' '3' '4'}, [-0.7 2.5], 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-200 0]);
    EEG.setname = [num2str(subi),'_PreICA'];
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), '_PreICA.set'],'filepath','F:\face-random\Preprocessing\PreICA');
end
%% Check the data and disgard bad trials.


%% Stage1-RunICA
clear all
clc
cd ('F:\face-random\Preprocessing\PreICA')

for subi=101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_PreICA.set'],'filepath','F:\face-random\Preprocessing\PreICA\');
    EEG = pop_runica(EEG,'icatype','fastica','approach','symm');
    EEG.setname = [num2str(subi),'_runICA'];
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), '_runICA.set'],'filepath','F:\face-random\Preprocessing\ICA');
end
% Stage1-ICA exclude (refer)
% 101 2; 102 2; 103 2; 104 7; 105 6; 106 8; 107 6; 108 1; 109 1; 110 6; 
% 111 1; 112 2; 113 1; 114 6; 115 2; 116 2; 117 1; 118 2; 119 9; 120 4; 121 2
%%
clear all
clc
cd ('F:\Face-random\Preprocessing\ICA')

for subi=101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_runICA.set'],'filepath','F:\face-random\Preprocessing\ICA\');
    EEG = pop_runica(EEG,'icatype','fastica','approach','symm');
    EEG.setname = [num2str(subi),'_runICA'];
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), '_runICA.set'],'filepath','F:\face-random\Preprocessing\ICA');
end
% sub 110 115 117 119 121 were excluded!!!!!!
%% Stage2-select the conditions
clear all
clc
cd ('F:\face-random\Preprocessing\ICA')

data(1).name = '_Neutral.set';data(2).name = '_Happy.set';data(3).name = '_N2H.set';data(4).name = '_H2N.set';
data(1).file = 'Condition1_Neutral'; data(2).file = 'Condition2_Happy';data(3).file = 'Condition3_N2H';data(4).file = 'Condition4_H2N';
for subi = 101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_runICA.set'],'filepath','F:\face-random\Preprocessing\ICA');
    EEG = pop_selectevent( EEG, 'type',1,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(1).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(1).name],'filepath',['F:\face-random\Preprocessing\Conditions\',data(1).file]);
end

for subi = 101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_runICA.set'],'filepath','F:\face-random\Preprocessing\ICA');
    EEG = pop_selectevent( EEG, 'type',2,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(2).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(2).name],'filepath',['F:\face-random\Preprocessing\Conditions\',data(2).file]);
end

for subi = 101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_runICA.set'],'filepath','F:\face-random\Preprocessing\ICA');
    EEG = pop_selectevent( EEG, 'type',3,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(3).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(3).name],'filepath',['F:\face-random\Preprocessing\Conditions\',data(3).file]);
end

for subi = 101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_runICA.set'],'filepath','F:\face-random\Preprocessing\ICA');
    EEG = pop_selectevent( EEG, 'type',4,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(4).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(4).name],'filepath',['F:\face-random\Preprocessing\Conditions\',data(4).file]);
end

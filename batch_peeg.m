function batch_lcnrl_eeg
% Batch file for the PEEG replication

mode = 'LU';

%%% Path %%%
if strcmp(mode,'LU')
    dir.dir = 'C:\Users\wurmf\Documents\Github\PEEG_tVNS\';
    %     dir.dir_data = 'C:\Users\wurmf\Documents\Data\PEEG_tVNS\'; %laptop
    dir.dir_data = 'C:\Users\wurmf\surfdrive\Shared\raw_data\EEG_data\'; %cloud
    chanFile = 'C:\Users\wurmf\Dropbox\MATLAB\Addons\ChanFiles\';
else
    warning('no mode specified! do you operate from ''home'' or are you at the ''LU''?')
end

%%% Directories %%%
% dir.BEHdir = strcat(dir.dir_data, 'data-BEH\');
dir.EEGdir = dir.dir_data; %strcat(dir.dir_data, 'data-EEG\');
dir.ERPdir = strcat(dir.dir_data, 'ERP\');
dir.TFAdir = strcat(dir.dir_data, 'TimeFreq\');
dir.mode = mode;

%%% Information about History %%%
% Change filename, but not structure name!
% For brevity, it will be referred to as INFO
fnhistory = 'PEEG_ExperimentInfo.mat';
dir.historyFile = fnhistory;
% alternative solution: global variables
% global fnhistory = 'PEEG_ExperimentInfo.mat'; %CAREFUL with global variables

%%% Add specific paths %%%
addpath('Routines\');
addpath('nway331-toolbox\')
addpath('eeg_2d_minimal_rejection\')

% Add fieltrip
ftpath = 'C:\Users\wurmf\Documents\Github\fieldtrip\';
% enableFieldTrip(1,ftpath) %TODO before using fieltrip routines
% addpath('C:\Users\wurmf\Documents\GitHub\fieldtrip\external\eeglab')
% enableFieldTrip(0,ftpath)

%%% Load participant structure %%%
% SCode = [1 2 3 4 6:18 20]; %full datasets
SCode = [4 6:18 20]; %good datasets
% % % SCode = [5 19]; % dropouts
% SCode = [8];

[S] = getS(SCode,dir);

%%% load electrode layout for this study (standard 32 electrode, T8 is replaced with FCz)
fid = fopen(fullfile(dir.dir,'elecLabels_64_standard.txt'));
elecinfo = textscan(fid,'%d %s','delimiter','\t');
fclose(fid);


% standard folder labels (wont apply to PEEG?)
%    e - after epoching
%    h - after electrode interpolation
%    c - after trial rejection
%    i - after ICA
%    r - after manual IC rejection
%    a - after automatic IC rejection

%% data import
% % eeg1_import_leiden(S,chanFile,'highpass',0.1,'in','\','relabel',elecinfo,'nochan',[71:79],'referenceElectrodes',[69 70]) %pilot only
% eeg1_import_leiden(S,chanFile,'highpass',0.1,'in','\','referenceElectrodes',[69 70])
% eeg1_import_leiden(S,chanFile,'highpass',0.1,'lowpass',40,'in','\','relabel',elecinfo)

% eeg1_import_leiden(S,chanFile,'in','\','referenceElectrodes',[69 70]) %quick

%% epoching
% eeg2_epoch_peeg(S,33,[-5 10],'out','stimLocked\e') %'baseline',[-200 0]
% eeg2_epoch_peeg(S,[45 46],[-1 4 1],'out','breakLocked\e') %'baseline',[-200 0]


%% preprocessing

%my standard method
% eeg3a_electrodeInterpolation(S,'kurtosis',5,'probability',5,'spectrum',3,'linkingType','or','in','stimLocked\e\','out','stimLocked\he\','electrodes',[1:64]);%[3:32 36:64]

% eeg3a_electrodeInterpolation(S,'kurtosis',5,'probability',5,'spectrum',3,'linkingType','or','in','breakLocked\e\','out','breakLocked\he\','electrodes',[1:64]);%[3:32 36:64]
% % % eeg3b_trialRejection(S,'absolute',300,'probability',5,'in','breakLocked\he\','out','breakLocked\che\','electrodes',[3:32 36:64]);

%sharons method (approximately)
% eeg3_preproc_sharon(S,'in','stimLocked\he','out','stimLocked\she','inspect','threshWindow',[0 3],'saveInfo','tvnsclean');
% eeg3_preproc_sharon(S,'in','breakLocked\he','out','breakLocked\she','threshWindow',[0 3],'reepoch',[0 3]);

%convert stimulation data to time-frequency
fft1_peeg(S,'stimLocked\she','stimLocked\she')

%extract baseline alpha from break data
% fft1_peeg_break(S,'breakLocked\she\','breakLocked\she\')
% fft2_peeg_parafac(S,'breakLocked\she\','CCD')
% fft2_peeg_parafac(S,'breakLocked\she\','PARAFAC')
% fft3_peeg_weight(S,'stimLocked\e','PARAFAC-results')
fft3_peeg_weight(S,'stimLocked\she','PARAFAC-results')




    function [S] = getS(SCode,dir)
        nS = length(SCode);
        for iVP = 1:nS
            
            S(iVP).index = SCode(iVP);
            S(iVP).dir = dir.dir;
            S(iVP).dir_data = dir.dir_data;
            S(iVP).EEGdir = dir.EEGdir;
            S(iVP).ERPdir = dir.ERPdir;
            S(iVP).TFAdir = dir.TFAdir;
            
            %             S(iVP).code = sprintf('Pilot%d',SCode(iVP));
            %             S(iVP).EEGfn = sprintf('Pilot%d',SCode(iVP));
            S(iVP).code = sprintf('P%03d',SCode(iVP));
            S(iVP).EEGfn = sprintf('P%03d',SCode(iVP));
            
            S(iVP).historyFile = dir.historyFile;
            
            S(iVP).suffix = '';
        end
        
    end




end
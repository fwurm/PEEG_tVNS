function fft1_peeg_break(S,inprefix,outprefix)

options.inprefix = inprefix;
options.outprefix = outprefix;

%% step 1: convert to frequency domain
for i = 1:length(S)
    subfun(S(i),options);
end



function subfun(S,options)

fprintf('=> Subject: %d\n', S.index);

%% import EEG
inpath = [S.EEGdir options.inprefix ''];
data = pop_loadset('filename',[S.EEGfn '.set'],'filepath',inpath);

%resample to 128 Hz
srate = 128;
data = pop_resample(data,srate);


%convert to doubles (from singles)
data.data = double(data.data);

%% convert to fieldtrip format
ftdata = eeglab2fieldtrip(data,'preprocessing','none');


%% time-frequency the data

freqs = [5:(1/3):15]; %target frequencies
times = data.times./1000; %target time window
winlength = 1; %taper window in seconds
% basewin = [-1 0]; %baseline window
% basetype = 'relchange'; %type of baseline

%from fieltrip online
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'EEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hamming';
cfg.pad          = 'nextpow2';
cfg.foi          = freqs;                   % analysis 5 to 15 Hz in steps of 0.33 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*winlength;  % length of time window = 3 sec
cfg.toi          = times;                   % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.keeptrials   = 'yes';
TFdata = ft_freqanalysis(cfg, ftdata);


%apply baseline
% cfg = [];
% cfg.baseline     = basewin;
% cfg.baselinetype = basetype;
% [bTFR_stim] = ft_freqbaseline(cfg, TFRdat);
% [bTFR_sham] = ft_freqbaseline(cfg, TFR_sham);



%save data
outpath = [S.TFAdir options.outprefix '\'];
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

% TFdata = bTFR_stim;
save(fullfile(outpath,[S.code '-break.mat']),'TFdata','-v7.3');
% TFdata = bTFR_sham;
% save(fullfile(outpath,[S.code '-sham.mat']),'TFdata','-v7.3');

%plot all electrodes
% cfg = [];
% cfg.showlabels   = 'yes';
% cfg.zlim         = [-1 2];
% ft_multiplotTFR(cfg, TFR_stim);
% ft_multiplotTFR(cfg, TFR_sham);

%plot single electrode
% channel = 'Pz';
% cfg = [];
% cfg.maskstyle    = 'saturation';
% cfg.zlim         = [-1 1];
% cfg.channel      = channel;
% ft_singleplotTFR(cfg, bTFR_stim);
% ft_singleplotTFR(cfg, bTFR_sham);

%plot differences in power
% cfg = [];
% cfg.maskstyle    = 'saturation';
% cfg.zlim         = [-0.5 0.5];
% cfg.channel      = channel;
% TFRdiff = bTFR_stim;
% TFRdiff.powspctrm = bTFR_stim.powspctrm - bTFR_sham.powspctrm;
% ft_singleplotTFR(cfg, TFRdiff);


function eeg3_preproc_sharon(S,varargin)


options.inprefix = [];
options.outprefix = [];
options.inspect = 0;
options.saveInfo = 0;
options.threshWindow = [0 1]; %default
options.reepoch = [];
options.tvnsclean = 0;

for i = 1:length(varargin)
    if strcmp(varargin{i}, 'in')
        if length(varargin)>i
            options.inprefix = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''in'' is not valid!');
        end
    end
    if strcmp(varargin{i}, 'out')
        if length(varargin)>i
            options.outprefix = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''out'' is not valid!');
        end
    end
    if strcmp(varargin{i}, 'threshWindow')
        if length(varargin)>i
            options.threshWindow = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''threshWindow'' is not valid!');
        end
    end
    if strcmp(varargin{i}, 'inspect')
        options.inspect = 1;
    end
    if strcmp(varargin{i}, 'saveInfo')
        options.saveInfo = 1;
    end
    if strcmp(varargin{i}, 'tvnsclean')
        options.tvnsclean = 1;
    end
    if strcmp(varargin{i}, 'reepoch')
        if length(varargin)>i
            options.reepoch = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''reepoch'' is not valid!');
        end
    end
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%      Step 3: Preprocessing      %');
disp('%          Sharon style           %');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');


for i = 1:length(S)
    subfun(S(i),options);
end

function subfun(S,options)

if options.saveInfo
    if exist(S.historyFile)
        load(S.historyFile); %load file to note down rejected trials (visual inspection of tvns artefact)
    else
        error('file not found')
    end
end
    

fprintf('=> Subject: %d\n', S.index);


%% import EEG
inpath = [S.EEGdir options.inprefix ''];
data = pop_loadset('filename',[S.EEGfn '.set'],'filepath',inpath);

%% get trial number (prior to preprocessing)
% load(S.historyFile); %load file to note down rejected trials (visual inspection of tvns artefact)
% INFO(S.index).trials_pre = data.trials;
% save(S.historyFile,'INFO')


%% linear detrend
for i = 1:data.trials 
    data.data(:,:,i) = detrend(data.data(:,:,i)')'; 
end

%% notch filter (line noise, 50Hz)

freqs = [50:50:100];
step1 = data;

for iF = 1:length(freqs)
    step1 = pop_eegfiltnew(step1,'locutoff',freqs(iF)-1,'hicutoff',freqs(iF)+1,'revfilt',1);
end

%chose example to plot
trial = 2;
elec = 17; %TP7

% figure
% 
% subplot(2,2,[1])
% plot(data.times,data.data(elec,:,trial))
% xlabel('Times [ms]')
% ylabel('Amplitude [\muV]')
% title('raw data')
% 
% X = fftshift(fft(data.data(elec,:,trial)));
% subplot(223)
% plot(f,abs(X)/N);
% xlabel('Frequency (in hertz)');
% title('No cleaning');
% 
% subplot(2,2,2)
% plot(step1.times,step1.data(elec,:,trial))
% xlabel('Times [ms]')
% ylabel('Amplitude [\muV]')
% title('raw data')
% 
% X = fftshift(fft(step1.data(elec,:,trial)));
% subplot(224)
% plot(f,abs(X)/N);
% xlabel('Frequency (in hertz)');
% title('pop-eegfiltnew notch [50:50:200]');


%% manual inspection
if options.inspect
    
    %first step: save epoched data (just to be save)
    outpath = [S.EEGdir options.outprefix '\'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end
    step1 = pop_saveset( step1, 'filename' , [S.EEGfn '.set'],'filepath',outpath);
    

    % step2: manual rejection using eegplot (load, reject, save)
    cmd = [ ...
        '[EEG] = pop_loadset(''filename'',''' [S.EEGfn '.set'] ''',''filepath'',''' outpath ''');'...
        '[tmprej tmprejE] = eegplot2trial( TMPREJ, ' num2str(step1.pnts) ', ' num2str(step1.trials) ');' ...
        'load(''' S.historyFile ''');' ...
        'INFO(' num2str(S.index) ').noArtefact = find(tmprej==1);' ...
        'save(''' S.historyFile ''',''INFO'');' ...
        '[EEGTMP LASTCOM] = pop_rejepoch(EEG, tmprej, 1);' ...
        '[EEG] = pop_saveset(EEGTMP,''filename'',''' [S.EEGfn '.set'] ''',''filepath'',''' outpath ''');' ...
        'clear EEGTMP;' ...
        ] ;
    eegplot(step1.data, ...
        'command',cmd, ...
        'butlabel','reject' ...
        );
    uiwait
    
    % step3: load the updated dataset
    step1 = pop_loadset('filename',[S.EEGfn '.set'],'filepath',outpath);
    
    if options.saveInfo
        INFO(S.index).nTrial_manualRef = step1.trials;
    end
else
    fprintf('skipping visual inspection...')
end


%% bandpass filter (5-15 Hz)
step2 = step1;

order = 3;
cutoff = [5 15];
nyquist = step2.srate/2;

[butterB,butterA] = butter(order,cutoff./nyquist,'bandpass'); %construct filter

% dat_orig = permute([data3.data],[2 1 3]); 
% dat_butt = filtfilt(butterB,butterA,double(dat_orig)); %apply filter
% data3.data = permute(single(dat_butt),[2 1 3]);

% figure

dat_orig = [step2.data];
dat_butt = nan(size(dat_orig));
for iChan = 1:size(dat_orig,1)
    dat_butt(iChan,:,:) = filtfilt(butterB,butterA,double(squeeze(dat_orig(iChan,:,:)))); %apply filter
%     plot(step2.times,dat_butt(iChan,:,2));
end
step2.data = single(dat_butt);

%% notch filter (tVNS noise, 25 Hz)
step3 = step2;
if options.tvnsclean
      
    freqs = [25:25:100];
    for iF = 1:length(freqs)
        step3 = pop_eegfiltnew(step3,'locutoff',freqs(iF)-1,'hicutoff',freqs(iF)+1,'revfilt',1);
    end
    
    % some plotting again
    trial = 2;
    elec = 17; %TP7
    
    Fs = data.srate;
    N = data.pnts;
    dF = Fs/data.pnts;
    f = -Fs/2:dF:Fs/2-dF;
    
%     figure
%     
%     subplot(2,3,[1 2 3])
%     hold on
%     diffsig = step3.data(elec,:,trial)-step2.data(elec,:,trial);
%     plot(step1.times,step1.data(elec,:,trial))
%     plot(step2.times,step2.data(elec,:,trial))
%     plot(step3.times,step3.data(elec,:,trial))
%     plot(step2.times,diffsig)
%     xlabel('Times [ms]')
%     ylabel('Amplitude [\muV]')
%     title('raw data')
%     legend({'original' 'bandpass' '25Hz+harm filt' 'diff band-25Hz'})
%     
%     X = fftshift(fft(step1.data(elec,:,trial)));
%     a = subplot(234)
%     plot(f,abs(X)/N);
%     xlabel('Frequency (in hertz)');
%     title('Before bandpass');
%     a.XLim = [0 250];
%     
%     X = fftshift(fft(step2.data(elec,:,trial)));
%     a = subplot(235)
%     plot(f,abs(X)/N);
%     xlabel('Frequency (in hertz)');
%     title('Before 25 Hz notch');
%     a.XLim = [0 250];
%     
%     X = fftshift(fft(step3.data(elec,:,trial)));
%     a = subplot(236)
%     plot(f,abs(X)/N);
%     xlabel('Frequency (in hertz)');
%     title('After 25 Hz notch');
%     a.XLim = [0 250];
else
    fprintf('skipping 25 Hz filtering (tvsn noise)\n')
end

%% reepoching (for break data)
if ~isempty(options.reepoch)
    step3 = pop_epoch(step3, {1}, options.reepoch,'epochinfo','yes'); 
    
    %delete duplicate events (end of each trial as they are overlapping)
%     [~,ia,~] = unique([step3.event.epoch],'first');
    a = {step3.epoch.event};
    b = cellfun(@length,a);
    c = ismember(b,1);
    step3 = pop_rejepoch(step3,~c,0);
end

%% thresholding

step4 = step3;

%method 1: standard EEGlab
thr = 200; 
chans = [1:64]; % default: without channels with blinks: 1 2 33 34 35
% chans2exclude = [7 8 15 16]; %close to left ear
% chans(ismember(chans,chans2exclude)) = [];

step4 = pop_eegthresh(step4,1,chans ,-thr,thr,options.threshWindow(1),options.threshWindow(2),0,1);
marks_th = step4.reject.rejthresh;

if options.saveInfo
    INFO(S.index).threshArtefact = find(marks_th==1);
    INFO(S.index).nTrial_thresh = step4.trials;
end

%method 2: Sharon
preEEGdata = [data.data];
tselect = dsearchn(data.times',options.threshWindow'.*1000);
EEGdata = preEEGdata(:,tselect(1):tselect(2),:);
cEEGdata = squeeze(mat2cell(EEGdata,size(EEGdata,1),size(EEGdata,2),ones(size(EEGdata,3),1)))';

% use a absloute crietria to reject trials/channels:
MAXABS_THRESHOLD=100;
max_abs_values=cell2mat(cellfun(@(x) max(abs(x),[],2),cEEGdata,'un',0)); % max_ abs values chan*trl
unvalid=max_abs_values>MAXABS_THRESHOLD; % this a boolean matrix of the currenly unvalid chan*trl combinations

view_on=0;
% [exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on,'trl');
[exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on,'chan');

%visualize sharons methods
figure
imagesc(unvalid)
ylabel('Channels')
xlabel('Trials')
b=get(gca,'XTickLabel');
set(gca,'XTickLabel',b,'fontSize',14)
for chan=exc_chan
text(1,chan,'*','Color','r','fontSize',30)
end
for trl=exc_trl
text(trl,1,'*','Color','r','fontSize',30)
end

%% interpolation?
warning('still need to setup interpolation routine')

%% save data
outpath = [S.EEGdir options.outprefix '\'];
if ~exist(outpath, 'dir')
    mkdir(outpath)
end
finaldata = pop_saveset( step4, 'filename' , [S.EEGfn '.set'],'filepath',outpath);

if options.saveInfo
    save(S.historyFile,'INFO')
end


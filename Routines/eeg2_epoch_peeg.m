function eeg2_epoch_peeg(S,targetEventI,epoch1,varargin)

warning on

options.epoch = epoch1;
options.baseline = [];
options.targetEvent = targetEventI;

options.inprefix = '';
options.outprefix = '';

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
    if strcmp(varargin{i}, 'baseline')
        if length(varargin)>i
            options.baseline = varargin{i+1};
            options.baseline_ind = 1;
        else
            disp('ERROR: Input for parameter ''baseline'' is not valid!');
        end
    end
end

%%%%%%%%%%%%%
% Portcodes %
%%%%%%%%%%%%%
trigger.stimon = [33];
trigger.stimoff = [34];
trigger.condtype1 = [31]; %condition types
trigger.condtype2 = [32];
trigger.breakon = [45];
trigger.breakoff = [46];
options.trigger = trigger;

if length(targetEventI)==1 && ~isempty(epoch1)
    options.epochtype = 'stimLocked';
elseif length(targetEventI)==2 && ~isempty(epoch1)
    options.epochtype = 'breakLocked';
else
    error('no correct epoching type can be specified\ncheck input parameters')
end
    




experiment.nBlock = 8;
experiment.allTrial = 88;
experiment.nTrial = 11;
% experiment.allRun = ;
% experiment.nRun = ;
options.experiment = experiment;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%      Step 2: epoching      %');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

fprintf('   Time-Locking Event: %s\n', options.targetEvent);
fprintf('   options.baseline: %s\n', mat2str(options.baseline));
fprintf('   options.epoch: %s\n', mat2str(options.epoch));


for i = 1:length(S)
    subfun(S(i),options);
end

function subfun(S,options)

cS = S.index; %current subject
fprintf('=> Subject: %d\n', cS);


%% get trial number (prior to preprocessing)
TEST.subject = cS;
TEST.expectTrial = nan;

%%%%%%%%%%%%%%%%%%
%expected events?%
%%%%%%%%%%%%%%%%%%
%Do the math! How many triggers would we expect if we assume perfect data
%(which unfortunately it wasnt...)
TEST.expectEvents= nan;


%% epoching

% import EEG
inpath = [S.EEGdir options.inprefix ''];
raw = pop_loadset('filename',[S.EEGfn '.set'],'filepath',inpath);

events = [raw.event];

if isstr(events(1).type)
    eventnums = cell2mat(cellfun(@(x) str2num(erase(x,'condition ')),{events.type},'UniformOutput',false));
else
    eventnums = [events.type];
end

TEST.imported = length(eventnums);



    
%% removing outliers (trigger offset due to external buttons)


% %option 1
% i = 0;
% prevents = eventnums;
% outliers = find(eventnums>255);
% while any(outliers) 
%     
%     i = i+1;
%     fprintf('   removing outliers round %d...\n',i)
%     
%     prevents(outliers) = prevents(outliers) - 8192;
%     
%     resolved = find(prevents==0);
%     events(resolved) = [];
%     prevents(resolved) = [];
%     
%     outliers = find(prevents>255);
% end
% eventnums = prevents;


%options 2
newenums = rem(eventnums,8192);
if ~isequal(rem(newenums,8192),eventnums)
    warning('new outlier detected')
end
eventnums = newenums;

TEST.outlierRem = length(eventnums);



% recoding of triggers (each experiment starts with 31 first)
trigger = options.trigger;
condtype = rem(S.index,2);
if (condtype == 1)
    trigger_tvns = trigger.condtype1;
    trigger_sham = trigger.condtype2;
elseif (condtype == 0)
    trigger_tvns = trigger.condtype2;
    trigger_sham = trigger.condtype1;
end
    

if (S.index == 9) %reversed conditions
    trigger_tvns = trigger.condtype2;
    trigger_sham = trigger.condtype1;
end



stimtrigger = [];
condtype = nan;
condcount = 0;
trialcount = 0;
warning('still need to fine-tune condition labelling during epoching')

if strcmp(options.epochtype,'stimLocked')
    for iE = 1:length(events)-1 
        
        %     if events(iE).type == 253
        %         condcount = condcount+1;
        %         condtype = rem(condcount,2); %1=tVNS, 0 = sham
        %     end  
        
        if eventnums(iE) == trigger_tvns
            condcount = condcount+1;
            condtype = 1; %tVNS,
        elseif eventnums(iE) == trigger_sham
            condcount = condcount+1;
            condtype = 0; %sham
        end        
        if eventnums(iE) == 33 %& eventnums(iE+1) == 34
            trialcount = trialcount+1;
            stimtrigger(trialcount).type = eventnums(iE);
            stimtrigger(trialcount).urevent = events(iE).urevent;
            stimtrigger(trialcount).latency =events(iE).latency;
            stimtrigger(trialcount).blocknum = condcount;
            stimtrigger(trialcount).trialnum = trialcount;
            stimtrigger(trialcount).condition = condtype;
        end
    end    
    TEST.expectBlock = nan;
    TEST.B_counted = condcount;       
    TEST.E_counted = length(stimtrigger);
       
    % identify events that are temporally (too) close
    timediff = diff([stimtrigger.latency]);
    tooshort = find(timediff<100); %sampling points?
    if ~isempty(tooshort)
        stimtrigger(tooshort+1) = []; %remove duplicate trials
    elseif (S.index == 16) %remove trials from short block
        stimtrigger([23 24 25]) = [];
    end
    
    
    
    raw.event = stimtrigger;
    TEST.E_duplicateRem = length(stimtrigger);
    
    
    if ~exist(S.historyFile)
        INFO = TEST;
    else
        load(S.historyFile); %load file to note down rejected trials (visual inspection of tvns artefact)
        INFO(S.index) = TEST;
    end    
    save(S.historyFile,'INFO')
    
elseif strcmp(options.epochtype,'breakLocked')
    
    
    
    gnu = 1;
    for iE = 1:length(events)-1        
        if eventnums(iE) == options.targetEvent(1)
            trialcount = trialcount+1;
            breakLatency(trialcount,1) = raw.event(iE).latency;
        end        
    end
    breakLatency(find(diff(breakLatency)<10*raw.srate))=[]; %clear duplicate events (closer than 10 seconds)
    
    breakends = find(eventnums==options.targetEvent(2));
    breakends(find(diff(breakends)<10))=[]; %clear duplicate events (more than 10 events inbetween)
    if length(breakends)==length(breakLatency)
        breakLatency(:,2) = [raw.event(breakends).latency];
    else
        error('break starts and ends still mismatch')
    end
    
    breaktrigger = [];
    breaktrial = 0;
    for iB = 1:size(breakLatency,1)
        currLat = breakLatency(iB,1);
        while (currLat < breakLatency(iB,2)-options.epoch(2)*raw.srate)
            breaktrial = breaktrial + 1;
            
            breaktrigger(breaktrial).type = 1;
            breaktrigger(breaktrial).urevent = breaktrial;
            breaktrigger(breaktrial).latency = currLat + (-1*options.epoch(1)*raw.srate); %
            breaktrigger(breaktrial).breaknum = iB;
            breaktrigger(breaktrial).trialnum = breaktrial;
            
            
            currLat = currLat+(sum(abs(options.epoch(1:2)))-options.epoch(3))*raw.srate;
        end
    end
    
    raw.event = breaktrigger;
    
%     EEG = eeg_regepochs(raw, 'recurrence',4, 'eventtype','3','extractepochs','off');
    
else
    error()
    
end



data = pop_epoch(raw,{},options.epoch(1:2), 'newname', [options.outprefix S.EEGfn '_epochs'], 'epochinfo', 'yes');

% if strcmp(options.epochtype,'breakLocked')
%     %delete duplicate events (end of each trial as they are overlapping)
%     [~,ia,~] = unique([data.event.urevent],'last');
%     data.event = data.event(ia);
% end


%% save data
outpath = [S.EEGdir options.outprefix '\'];
if ~exist(outpath, 'dir')
    mkdir(outpath)
end
finaldata = pop_saveset( data, 'filename' , [S.EEGfn '.set'],'filepath',outpath);

% load('rejectedTrials.mat'); %load file to note down rejected trials (visual inspection of tvns artefact)
% rejectedTrials(S.index).trials_pre = length(stimtrigger);
% save('rejectedTrials.mat','rejectedTrials')




























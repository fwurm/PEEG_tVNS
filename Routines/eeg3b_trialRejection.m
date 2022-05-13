function eeg3b_trialRejection(S,varargin)
% eeg3b_trialRejection() - Removes epochs (trials) with abnormal data.
%
% Usage: 
%   >> eeg3b_trialRejection(S, optionalParameters)
% 
%   Mandatory Inputs:
%           S           - Information on the subject sample (a structure
%                         that has to be created in the batch file)              
%
%   Optional Inputs:
%           'absolute'  - Removes electrodes based on an absolute threshold
%                         criterium (e.g., 300 ?V).
%           'probability' - Removes electrodes based on the joint
%                         probability criterium. Specify the number of 
%                         standard deviations from the mean. (e.g., 5)
%           'linkingType' - Specify whether an electrode has to be strange 
%                         according to all above selected criteria ('and')
%                         or any one of the criteria ('or'). Default: 'or'
%           'electrodes'- All electrodes that are considered for
%                         interpolation. (Default: [3:32 36:64], i.e., all
%                         electrodes of the BioSemi ActiveTwo system except
%                         the first row, because eye blinks will normally 
%                         be corrected via ICA)
%           'blinks'    - Specify a time interval from which trials with
%                         eye blinks are to be removed (e.g. [0 200 500]
%                         to reject blinks from 0-200 ms with an offset
%                         period of 500ms) BETA! MIGHT NOT WORK PROPERLY!
%           'artefacts' - Specify the name of a file that contains all 
%                         information on artefact rejection and correction
%                         (Default: ''ArtefactCorrectionData.mat').
%           'in'        - folder name and filename prefix for input eeg
%                         files to be loaded. Low dashes ("_") indicate 
%                         new folders ("\")
%           'out'       - folder name and filename prefix for output eeg
%                         files to be saved. Low dashes ("_") indicate 
%                         new folders ("\")
%
%   Example: 
%   >> eeg3b_trialRejection(S(cS), 'threshold', 300, 'probability', 5,...
%           'linkingType', 'or', 'electrodes', [1:64], 'artefacts', 'MyArtefacts.mat',...
%           'in', 'Filtered01-40_c1e', 'out', 'Filtered01-40_c2e' );
%                                         
% (c) 2017 by Robert Steinhauser



inprefix = '';
outprefix = '';
blinkRejectionInterval = [0 0 0];
horzEOGInterval = [0 0 0];
noPlot = 0;

absoluteThreshold = [];
probabilityThreshold = [];
channelsWithoutBlinks = [3:32 36:64]; % default: without channels with blinks: 1 2 33 34 35
artefactDataFile = 'ArtefactCorrectionData.mat';
rejWindow = [];


for i = 1:length(varargin)
    if strcmp(varargin{i}, 'in') 
       if length(varargin)>i  
           inprefix = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''in'' is not valid!');
       end
    end   
    if strcmp(varargin{i}, 'out') 
       if length(varargin)>i  
           outprefix = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''out'' is not valid!');
       end
    end        
    if strcmp(varargin{i}, 'probability') 
       if length(varargin)>i && isnumeric(varargin{i+1})  
           probabilityThreshold = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''probabilityThreshold'' is not valid!');
       end
    end          
    if strcmp(varargin{i}, 'absolute') 
       if length(varargin)>i && isnumeric(varargin{i+1})  
           absoluteThreshold = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''absoluteThreshold'' is not valid!');
       end
    end    
    if strcmp(varargin{i}, 'blinks') 
       if length(varargin)>i  
           blinkRejectionInterval = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''blinks'' is not valid!');
       end
    end    
    if strcmp(varargin{i}, 'linkingType')
        if length(varargin)>i && (strcmp(varargin{i+1},'and') || strcmp(varargin{i+1},'or'))
            linkingType = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''linkingType'' is not valid!');
        end
    end
    if strcmp(varargin{i}, 'electrodes')
        if length(varargin)>i && isnumeric(varargin{i+1}) 
            channelsWithoutBlinks = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''electrodes'' is not valid!');
        end
    end   
    if strcmp(varargin{i}, 'artefacts')
        if length(varargin)>i
            artefactDataFile = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''artefacts'' is not valid!');
        end
    end  
    if strcmp(varargin{i}, 'rejWindow')
        if length(varargin)>i
            rejWindow = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''rejWindow'' is not valid!');
        end
    end
    

end





%%%%%%%%%
% START %
%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%          Step 3b: Trial Rejection          %');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');




fprintf('   Subject\tTotTrials\tTotBlinks\tTotHorzEOG\tBlinks(%s-%s)\tHorzEOG(%s-%s)\tOutliers\tTrialsWithoutBlinks\tTrialsWithoutOutliers',...
    blinkRejectionInterval(1) ,blinkRejectionInterval(2) , horzEOGInterval(1),horzEOGInterval(2));

% remove artifacts
for i = 1:length(S)
    out = subfun2(S(i),i,inprefix,outprefix, probabilityThreshold, absoluteThreshold, channelsWithoutBlinks, blinkRejectionInterval,artefactDataFile,rejWindow);
   % for j = 1:length(out)
    %    outall{j}(i) = out{j};
   % end
end




function out = subfun2(S,ind,inprefix,outprefix, probabilityThreshold, absoluteThreshold, channelsWithoutBlinks, blinkRejectionInterval,artefactDataFile,rejWindow)


 

inpath = [S.EEGdir inprefix '\'];

data = pop_loadset('filename', [S.EEGfn '.set'], 'filepath', inpath);

rejTrials = zeros(size(data.data,3),1);

trialsTot = size(data.data,3);

blinkCountTotal = 0;
blinkCountTimeWindow = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%     BLINK rejection    %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if blinkRejectionInterval ~= [0 0 0]

    rawBlinks = data;
    
    blinkData = squeeze(rawBlinks.data(67,:)-rawBlinks.data(68,:));
    
    blinkData(blinkData<150)=0;
    blinkData(blinkData>=150)=200;
    
    rawBlinks.data(73,:) = blinkData;
    rawBlinks.nbchan = size(rawBlinks.data,1);
    rawBlinks.chanlocs(73).labels = 'blinks';
    
    blinkCountTotal = 0;
    blinkCountTimeWindow = 0;
    
    
    for i=1:size(data.data,3)
        
        m = squeeze(rawBlinks.data(73,:,i));
        if numel(m(m>0))>0
            %disp(['Trial mit Blink: ' num2str(i)]);
            blinkCountTotal = blinkCountTotal+1;
        end
        
        
        timeStart = rs_TimeToFrame(blinkRejectionInterval(1),512,blinkRejectionInterval(3)); % 800; % ms
        timeEnd = rs_TimeToFrame(blinkRejectionInterval(2),512,blinkRejectionInterval(3)); % 1250; % ms
        
        m = squeeze(rawBlinks.data(73,timeStart:timeEnd,i));
        if numel(m(m>0))>0
            blinkCountTimeWindow = blinkCountTimeWindow+1;
            rejTrials(i) = 1;
        end
        
    end
    
    
    
    
data = pop_rejepoch(data, rejTrials,0);

end
 %%% END BLINK-STUFF
 

 
 trialsAfter1 = size(data.data,3);
  
  
 

 
 
 
%% remove artifacts


chan = channelsWithoutBlinks; % without channels with blinks: 1 2 33 34 35
marks_th = [];
% reject epochs with outlier values over thr mV 
if ~isempty(absoluteThreshold)
    thr = absoluteThreshold;
    
    starttime = data.xmin;
    endtime = data.xmax;
    if ~isempty(rejWindow)
        starttime = rejWindow(1);
        endtime = rejWindow(2);
    end
    
    data = pop_eegthresh(data,1,[chan] ,-thr,thr,starttime,endtime,0,0); %,-1,0.998,2,1);
    marks_th = data.reject.rejthresh;
end

% reject epochs with improbable data
marks_jp = [];
if ~isempty(probabilityThreshold)
%     try
        data = pop_jointprob(data,1,[chan] ,probabilityThreshold,probabilityThreshold,0,0);
        
%     catch
%         disp('ERROR');
%     end
    
    marks_jp = data.reject.rejjp;
    
    marks_jp = marks_jp*2;
else
    marks_jp = zeros(1,size(data.data,3));    
end

marks = marks_th + marks_jp;


% rt_disp(rt_sortbynum(rtcell(rt_condfreq(marks)),1))

% save([S.EEGfn '_rejep.mat'], 'marks');

out{1} = {S.EEGfn};
out{2} = length(find(marks==1));
out{3} = length(find(marks==2));
out{4} = length(find(marks==3));

marks(marks>0) = 1;
data2 = data;
data = pop_rejepoch(data, marks, 0);
outliers = length(find(marks));

EOGCountTotal = 0;
EOGCountTimeWindow = 0;


fprintf('   %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\n',...
    S.index,trialsTot,blinkCountTotal,EOGCountTotal,blinkCountTimeWindow,EOGCountTimeWindow,outliers,trialsAfter1,size(data.data,3))


 warning = 0;
 if outliers > 100    
     pop_eegplot(data2,1,0,0);
     fprintf('WARNING: More than 100 Outliers with VP %d\n', S.index);
     warning = 1;
 end
 
 
 
 
 
 

if exist(artefactDataFile, 'file') == 2
    load(artefactDataFile);
else
    artefactData = [];
end

artefactData.TrialRejection{S.condnum}(S.index).outliersTogether = outliers;
artefactData.TrialRejection{S.condnum}(S.index).WARNING = warning;

artefactData.TrialRejection{S.condnum}(S.index).byAbsoluteThreshold = length(find(marks==1));
artefactData.TrialRejection{S.condnum}(S.index).byProbability = length(find(marks==2));
artefactData.TrialRejection{S.condnum}(S.index).byBlinks = blinkCountTotal;


artefactData.TrialRejection{S.condnum}(S.index).trialsBefore = trialsTot;
artefactData.TrialRejection{S.condnum}(S.index).trialsAfterwards = size(data.data,3);
artefactData.TrialRejection{S.condnum}(S.index).rejectedTrials = {marks};


save(artefactDataFile,'artefactData');
 
 
outpath = [S.EEGdir outprefix '\'];
if ~exist(outpath, 'dir')
       mkdir(outpath)
end
 
pop_saveset(data,'filename', [S.EEGfn '.set'], 'filepath', outpath);





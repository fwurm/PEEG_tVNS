function eeg3a_electrodeInterpolation(S,varargin)
% eeg3a_electrodeInterpolation() - Removes electrodes with abnormal data and 
%                 subsequently conducts a spherical interpolation with 
%                 information from the surrounding electrodes. 
%
% Usage: 
%   >> eeg3a_electrodeInterpolation(S, optionalParameters)
% 
%   Mandatory Inputs:
%           S           - Information on the subject sample (a structure
%                         that has to be created in the batch file)              
%
%   Optional Inputs:
%           'kurtosis'  - Removes electrodes based on the kurtosis
%                         criterium. Specify the number of standard
%                         deviations from the mean. (e.g., 5)
%           'probability' - Removes electrodes based on the joint
%                         probability criterium. Specify the number of 
%                         standard deviations from the mean. (e.g., 5)
%           'spectrum'  - Removes electrodes based on the spectrum
%                         criterium. Specify the number of whatever
%           'linkingType' - Specify whether an electrode has to be strange 
%                         according to all above selected criteria ('and')
%                         or any one of the criteria ('or'). Default: 'or'
%           'electrodes'- All electrodes that are considered for
%                         interpolation. (Default: [3:32 36:64], i.e., all
%                         electrodes of the BioSemi ActiveTwo system except
%                         the first row, because eye blinks will normally 
%                         be corrected via ICA)
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
%   >> eeg3a_electrodeInterpolation(S(cS), 'kurtosis', 5, 'probability', 5,...
%           'linkingType', 'or', 'electrodes', [1:64], 'artefacts', 'MyArtefacts.mat',...
%           'in', 'Filtered01-40_e', 'out', 'Filtered01-40_c1e' );
%                                         
% (c) 2017 by Robert Steinhauser
 
 

inprefix = '';
outprefix = '';



kurtosisThreshold = [];
probabilityThreshold = [];
spectrumThreshold = [];
linkingType = 'and';
channelsWithoutBlinks = [3:32 36:64]; % default: without channels with blinks: 1 2 33 34 35
artefactDataFile = 'ArtefactCorrectionData.mat';

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
    if strcmp(varargin{i}, 'kurtosis') 
       if length(varargin)>i && isnumeric(varargin{i+1})  
           kurtosisThreshold = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''kurtosisThreshold'' is not valid!');
       end
    end      
    if strcmp(varargin{i}, 'probability') 
       if length(varargin)>i && isnumeric(varargin{i+1})  
           probabilityThreshold = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''probabilityThreshold'' is not valid!');
       end
    end          
    if strcmp(varargin{i}, 'spectrum') 
       if length(varargin)>i && isnumeric(varargin{i+1})  
           spectrumThreshold = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''spectrumThreshold'' is not valid!');
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
            disp('ERROR: Input for parameter ''artefactDataFile'' is not valid!');
        end
    end         
     
    
end





%%%%%%%%%
% START %
%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%      Step 3a: Electrode Interpolation      %');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');


% interpolate channels
for i = 1:length(S)
    if (S(i).index) == 10 && (S(i).session ==2)
        gnu = 1;
    end
    out = subfun1(S(i),i,inprefix,outprefix,kurtosisThreshold, probabilityThreshold, spectrumThreshold, linkingType, channelsWithoutBlinks,artefactDataFile);
end




function out = subfun1(S,ind,inprefix,outprefix,kurtosisThreshold, probabilityThreshold, spectrumThreshold, linkingType, channelsWithoutBlinks,artefactDataFile)



 
 
 
%% interpolate channels

disp(['inprefix = ' inprefix]);

inpath = [S.EEGdir inprefix '\'];
data = pop_loadset( 'filename', [S.EEGfn '.set'], 'filepath', inpath);


dataOld{ind} = data;

chan = channelsWithoutBlinks; 


indelecp = [];
indeleck = [];
indelecs = [];

if ~isempty(probabilityThreshold)

    % probability
    [~, indelecp, measure] = pop_rejchan(data, 'elec', chan, 'threshold', probabilityThreshold, 'measure', 'prob', 'norm', 'on');
    indelecp = chan(indelecp); % !!!! bug fix !!!!
    e = [];
    for i = 1:length(indelecp)
        e = [e ' ' data.chanlocs(indelecp(i)).labels];
    end
    
    % if noPlot == 0;
    %     subplot(1,2,1); plot(sort(measure),'X');
    %     if ~isempty(e)
    %         text(10,2, e);
    %     end
    % end
    
    
end

if ~isempty(kurtosisThreshold)
    
    % kurtosis
    [~, indeleck, measure] = pop_rejchan(data, 'elec', chan, 'threshold', kurtosisThreshold, 'measure', 'kurt', 'norm', 'on');
    indeleck = chan(indeleck); % !!!! bug fix !!!!
    e = [];
    for i = 1:length(indeleck)
        e = [e ' ' data.chanlocs(indeleck(i)).labels];
    end
    
    % if noPlot == 0;
    %     subplot(1,2,2); plot(sort(measure),'X');
    %     if ~isempty(e)
    %         text(10,2, e);
    %     end
    % end;
    
end


if ~isempty(spectrumThreshold)

    %Spectrum
    [~, indelecs] = pop_rejchan(data,'elec',chan,'threshold',spectrumThreshold,'norm','on','measure','spec');
    indelecs = chan(indelecs); % !!!! bug fix !!!!
    
    e = [];
    for i = 1:length(indelecs)
        e = [e ' ' data.chanlocs(indelecs(i)).labels];
    end
    
    % if noPlot == 0;
    %     subplot(1,2,2);% plot(sort(measure),'X');
    %     if ~isempty(e)
    %         text(10,2, e);
    %     end
    % end;

end


% rejections

if strcmp(linkingType,'or')
    
    rc = [];
    for i = 1:64
        if (~isempty(find(indelecp==i, 1))) || ~isempty(find(indeleck==i, 1)) || ~isempty(find(indelecs==i, 1))
            rc = [rc i];
            
        end
    end
    
elseif strcmp(linkingType,'and')
    rc = [];
    if kurtosisThreshold == [],    indeleck = [1:64]; end
    if probabilityThreshold == [], indelecp = [1:64]; end
    if spectrumThreshold == [],    indelecs = [1:64]; end
    
    for i = 1:64
        if (~isempty(find(indelecp==i, 1))) && ~isempty(find(indeleck==i, 1)) && ~isempty(find(indelecs==i, 1))
            rc = [rc i];
            
        end
    end
    
else
    error('No valid linkingType in ElectrodeInterpolation!');
end

%data.badchannels = rc;


e = [];
for i = 1:length(rc)
    e = [e ' ' data.chanlocs(rc(i)).labels];
end
%set(gcf,'Name',[S.EEGfn ' - ' e]);



out{1} = {S.EEGfn};
out{2} = length(rc);

% save([S.EEGfn '_rejch.mat'],'rc');

s=[];
for i=1:length(rc)
    s = [s ' ' data.chanlocs(rc(i)).labels ];
end

fprintf('   Subject: %d - Channels: %s\n', S.index, s);

% apply interpolation


if ~isempty(indelecs)
    gnu = 1;
end
       


%pop_eegplot(data, 1, 0, 0);

if ~isempty(rc)
    data = pop_interp(data,rc,'spherical');
end

%pop_eegplot(data, 1, 0, 0);




if exist(artefactDataFile, 'file') == 2
    load(artefactDataFile);
else
    artefactData = [];
end

if length(rc)>0
    artefactData.ElectrodeInterpolation.Together{S.index,S.condnum} = data.chanlocs(rc).labels ;
else
    artefactData.ElectrodeInterpolation.Together{S.index,S.condnum} = [] ;
end
artefactData.ElectrodeInterpolation.Together{S.index,S.condnum} = rc;
artefactData.ElectrodeInterpolation.byProbability{S.index,S.condnum} = indelecp;
artefactData.ElectrodeInterpolation.byKurstosis{S.index,S.condnum} = indeleck;
artefactData.ElectrodeInterpolation.bySpectrum{S.index,S.condnum} = indelecs;


save(artefactDataFile,'artefactData');


outpath = [S.EEGdir outprefix '\'];
if ~exist(outpath, 'dir')
       mkdir(outpath)
end

pop_saveset(data,'filename', [S.EEGfn '.set'], 'filepath', outpath);





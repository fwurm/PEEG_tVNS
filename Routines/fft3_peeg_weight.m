function fft3_peeg_weight(S,indata,inparafac)

basewin = [-1 0]; %basline window as specified in Sharon

%components which reflect the baseline alpha topographie and spectrum
nAlpha = [4,1; 6,1; 7,2; 8,2; 9,1; 10,1; 11,1; 12,1; 13,1; 14,1; 15,1; 16,1; 17,1; 18,1; 20,1]; 

allfacs = [];
alldata = [];

%%loading relevant info
for i = 1:length(S)
    
    fprintf('   loading %s...\n',S(i).code)
    
    fin = fullfile(inparafac,[S(i).code '_PARAFAC.mat']);
    load(fin); %loaded as 'Factors' (dims = [trials chans freqs timepnts])
    allfacs{i,1} = Factors;
    
    if ~isempty(find(S(i).index==nAlpha(:,1)))
        allfacs{i,2} = nAlpha(find(nAlpha(:,1)==S(i).index),2);
    else
        error('participant is not parafact yet!')
    end
    
    
    outpath = [S(i).TFAdir indata '\'];
    load(fullfile(outpath,[S(i).code '-stim.mat'])); %loaded as 'TFdata'
    alldata{i,1} = TFdata;
    load(fullfile(outpath,[S(i).code '-sham.mat'])); %loaded as 'TFdata'
    alldata{i,2} = TFdata;
end

clear Factors
clear TFdata

%% applying PARAFAC weighting for single-channel (figure 4D in Sharon)
for i = 1:length(S)
    
    fprintf('   processing %s...\n',S(i).code)
    
    chanweight = allfacs{i,1}{2}(:,allfacs{i,2});
    freqweight = allfacs{i,1}{3}(:,allfacs{i,2});
    
    for iC = 1:2 %condition
        
        data = [alldata{i,iC}.powspctrm]; %4-D (trial,channel,frequency,time)
        
        sdata = size(data);
        
        %multiply with channelweights
        D = permute(data,[2 3 1 4]);
        E = reshape(D,sdata(2),prod(sdata([1 3 4])));
        F = chanweight'*E;
        data2 = reshape(F,sdata([1 3 4])); %3-D (trial,frequency,time)
        
        %multiply with freqweights
        D = permute(data2,[2 1 3]);
        E = reshape(D,sdata(3),prod(sdata([1 4])));
        F = freqweight'*E;
        data3 = reshape(F,sdata([1 4])); %2-D (trial,time)
        
        %baselining
        times = alldata{i,1}.time;
        bpts = dsearchn(times',basewin');
        basepwr = mean(data3(:,bpts(1):bpts(2)),2);
        
        data_bl = data3 ./ repmat(basepwr,[1 sdata(4)]);
        mdata{iC}(i,:) = mean(data_bl); %1-D (time)
    end
    
end

figure
a = axes;
plot(times,mean(mdata{1}),times,mean(mdata{2}))
a.XLim = [-1 6];
a.YLim = [0.7 1.7];
legend({'stim' 'sham'})


%plot with confidence interval (TODO: check validity of CIs)
%
% normd = [];
% ci = [];
% for j = 1:size(mdata{iC},2) % sample points
%     dall = (nanmean(mdata{1}(:,j))+nanmean(mdata{2}(:,j)))/2;
%     for iC = 1:2 % cons
%         for k = 1:size(mdata{iC},1) % subjects
%             meand(iC,:) = mean(mdata{iC},1);
%             dk = (mean(mdata{1}(k,j)) + mean(mdata{2}(k,j)))/2;
%             normd(iC,j,k) = mdata{iC}(k,j) - dk + dall;
%         end
%         ci(iC,j) = 1.96 * std(normd(iC,j,:))/sqrt(size(mdata{iC},2));
%     end
% end
% 
% figure;
% a = axes;
% hold on;
% bh = shadedErrorBar(times,meand(1,:),ci(1,:),'b',1);
% rh = shadedErrorBar(times,meand(2,:),ci(2,:),'r',1);
% a.XLim = [-1 6];
% a.YLim = [0.7 1.7];
% l = legend([bh.mainLine rh.mainLine],{'stim' 'sham'});


%% grandaverage (without any application of channel or frequency weighting)
pregrandavg = [];
cfg = [];
evalString = cell(1,2);
for i = 1:length(S)
    for iC = 1:2
        pregrandavg{i,iC} = ft_freqdescriptives(cfg,alldata{i,iC} );
        evalString{iC} = [evalString{iC} 'pregrandavg{' num2str(i) ',' num2str(iC) '},'];
    end
end
evalString{1}(end) = [];
evalString{2}(end) = [];



cfg.keepindividual = 'yes'; % 'yes' or 'no' (default = 'no')   HIER WAR "yes" von Clara ???
cfg.foilim         = 'all'; % [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
cfg.toilim         = 'all'; % [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
cfg.channel        = 'all'; % Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
cfg.parameter      = 'powspctrm'; %string or cell-array of strings indicating which parameter(s) to average. default is set to 'powspctrm', if it is present in the data.
    
eval(['[allgrandavg{1}] = ft_freqgrandaverage(cfg, ' evalString{1} ');']);
eval(['[allgrandavg{2}] = ft_freqgrandaverage(cfg, ' evalString{2} ');']);


clear alldata
clear pregrandavg
clear data data2 data3 D E F S





%% permmutation testing (without any application of channel or frequency weighting)
cfg = [];
cfg.channel          = 'EEG';
cfg.latency          = [0 5];
cfg.frequency        = 'all';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; % 'ft_statfun_actvsblT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum'; %%' maxsum', 'maxsize', 'wcm'
cfg.minnbchan        = 3; 


cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

%load channel configuration
elecFileName = [S(1).EEGdir '\stimLocked\e\' S(1).EEGfn '.set']; 
cfg_neighb.elecfile = elecFileName;
cfg_neighb.method    = 'triangulation';
data = ft_read_sens(elecFileName);
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, data);

%construct design matrix
subj = size(allgrandavg{1}.powspctrm,1);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

% perform cluster-based permutation
[stat] = ft_freqstatistics(cfg,allgrandavg{1},allgrandavg{2});

gnu = 1;

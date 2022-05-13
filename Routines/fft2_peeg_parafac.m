function fft2_peeg_parafac(S,inprefix,mode)

options.inprefix = inprefix;
options.mode = [];

if strcmp(mode,'CCD') %estimate number of components via core consistency diagnostic
    options.mode = 1;   
elseif strcmp(mode,'PARAFAC') % use estimated number to calculate PARAFAC model
    options.mode = 2;
else
    error('no mode specified')
end

%PARAFAC components, visually derived from core consistency plot
nComp = [4,2; 6,1; 7,3; 8,3; 9,3; 10,1; 11,2; 12,2; 13,2; 14,3; 15,2; 16,2; 17,2; 18,2; 20,2];
%in Sharon: "The proper number of components was determined by using the CCD, in which the number of components is the highest when the minimal value of CCD is 55%"
%in nway documentation: model is valid with 40%+




%% step 1: convert to frequency domain
for i = 1:length(S)
    
    % get number of components for the PARAFAC model
    if ~isempty(find(S(i).index==nComp(:,1)))
        sComp = nComp(find(nComp(:,1)==S(i).index),2);
    else
        warning('number of components is not determined yet')
    end
    
    subfun(S(i),sComp,options);
end

function subfun(S,nFactors,options)

close all

fprintf('\n==> subject %d\n',S.index)

%load time-frequency data
inpath = [S.TFAdir options.inprefix '\'];
fin = fullfile(inpath,[S.code '-break.mat']);
load(fin) %as TFdata

X = TFdata.powspctrm;

if options.mode == 1 %determine number of components from CCD criterion
    
    [ssX,Corco] = pftest(3,X,5,[0 0 0 0 NaN]);
    
    fout = strcat('PARAFAC-results\',S.code,'_CCD.fig');
    savefig(gcf,fout); %save figure
    fout = strcat('PARAFAC-results\',S.code,'_CCD.mat');
    save(fout,'ssX','Corco'); %save data from CCD test
    
elseif options.mode == 2 %fit PARAFAC model using the number of components determines by the CCD criterion
    
    opts = [1e-6 1 0 0 10 2500]; %options
    const = [2 2 2 2]; %constraints
    Factors = parafac(X,nFactors,opts,const);
    
    %load EEGLAB data for plotting only
    inpath = [S.EEGdir options.inprefix ''];
    data = pop_loadset('filename',[S.EEGfn '.set'],'filepath',inpath);    
    
    
    % plot results from PARAFAC (spatial and frequency domain)
    f = figure;
    for iF = 1:nFactors
        maxLoad = max(Factors{1,2}(:,iF));
        subplot(2,nFactors,iF)
        topoplotMS(Factors{1,2}(:,iF), data.chanlocs(1:64),'maplimits',[0 maxLoad]);
        title(sprintf('Component %d',iF))
    end
    subplot(2,nFactors,[nFactors+1:2*nFactors])
    plot(TFdata.freq,Factors{3})
    lgdtxt = cellfun(@(x) ['Component ' num2str(x)],num2cell([1:nFactors]),'UniformOutput',false);
    legend(lgdtxt)
    title('Frequency spectrum')
        
    fout = strcat('PARAFAC-results\',S.code,'_PARAFAC.fig');
    savefig(f,fout); %save figure
    fout = strcat('PARAFAC-results\',S.code,'_PARAFAC.mat');
    save(fout,'Factors'); %save data from PARAFAC
    
else
    error('no mode specified for parafac script')
end

gnu = 1;


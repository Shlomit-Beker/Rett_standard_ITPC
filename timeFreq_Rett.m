%% calculate TFR and ITC for one channel a time, averaged across subjects. 
%Used for Tufi's rett data
% RUN AFTER PL_Rett to get AllData structure, or upload it straight from here 
% Shlomit Beker shlomitbeker@gmail.com

% calculate ITC

clear dataTemp spectAll spectrumEst itc Group data DATA_I  itc_1 itcAllRett itcAllTd itcDiff 
addpath('C:\Users\sbeker\Dropbox (EinsteinMed)\GENERAL ANALYSIS FILES');
addpath('/Users/shlomit/Dropbox (EinsteinMed)/GENERAL ANALYSIS FILES');
addpath('/Users/shlomit/Dropbox (EinsteinMed)/TUFI/Rett_ITPC');

if  ~exist('AllData')
    load AllData
    AllData_ds = All_Data_ds(AllData); % calling function to down sample the data. 
else
    AllData_ds = All_Data_ds(AllData);
end

%%
Fs = 256;
chnI = 1;
gwidthWlt = 3;
freqoi = 2:25;
%freqoi = logspace(0.2,1.5,20);
widthWlt = linspace(3,5,30);
chns = 47;
%numTrl = 1; % enter the number of trials
t = -1:1/Fs:1;
 
% calculate ITC
% try one to initialize size of timeoi and freqoi
% data is a fieldtrip structure with single trial data
% inputs: 
% dataTemp - single trial data matrix (chans X time)
% t - time vector in seconds

%TD: 1-24; RETT: 25-42; % cohort in Tufi's Rett paper. 
DATA_I{1} = AllData_ds(:,1:24);
DATA_I{2} = AllData_ds(:,25:end);

%%
trialNumber_all(:,1)= sum(trialNumber(:,1:24),2);
trialNumber_all(:,2)= sum(trialNumber(:,25:end),2);

%%
Group = 1; %Change before run: 1-TD/2-RETT
if Group == 2
prompt = 'Remove noisy subjects? (1-Yes,0-No) ';
Remove = input(prompt);

if Remove == 1

    DATA_I{2}(:,1) = []; %bad subject#1
    DATA_I{2}(:,11) = [];%bad subject#12
    trialNumber_all(:,2)= trialNumber_all(:,2)-trialNumber(:,25)-trialNumber(:,36); %bad subject#1
end
    
end  
    
tic
dataTemp = DATA_I{Group}{1,1}(chns,:,1);
[~,freqoi,timeoi] = ft_specest_wavelet(dataTemp, t, 'freqoi', freqoi,'width', widthWlt, 'gwidth',gwidthWlt);

 for k =  1:size(DATA_I{Group},1) %conditions
     for m = 1:length(DATA_I{Group}(k,:)); % participants
        % initialize matrix for ITC
        % initialize matrix for spectrogram
        % this matrix holds the instantaneous phase for each trial, frequency
        % and time point 
        numTrl = size(DATA_I{Group}{k,m},3);
        spectAll = zeros(numTrl,length(freqoi),length(timeoi));
        for trlI = 1 : numTrl
            % select the correct trial and channel
            %dataTemp = DATA_I{Group}{k,m}(chns(chnI),:,trlI);
            dataTemp = DATA_I{Group}{k,m}(chns(chnI),:,trlI);

            % calculate time frequency analysis using wavelets
            [spectrumEst,freqoi,timeoi] = ft_specest_wavelet(dataTemp, t, 'freqoi', freqoi, 'width', widthWlt, 'gwidth',gwidthWlt);
            spectAll(trlI,:,:) = spectrumEst(1,:,:);
        end
    % calculate itc for a single channel across trials
    itc{m}(chnI,:,:) = it_calcITC(spectAll);
     end
   itc_1{k} = itc;

 end
 toc
 
 for k = 1:length(itc_1)
     sumItc = zeros(length(freqoi),length(t));
     for s = 1:length(itc_1{k})
         sqItc = squeeze(itc_1{k}{s});
         sumItc = sumItc+sqItc;
     end    
         if Group == 1
             itcAllTd{k} = sumItc./length(DATA_I{Group});
             data{k} = itcAllTd{k};
             itcTD{k} = cell2mat(itc_1{k}'); %concatenated subjects
         elseif Group == 2
             itcAllRett{k} = sumItc./length(DATA_I{Group});
             data{k} = itcAllRett{k};
             itcRETT{k} = cell2mat(itc_1{k}');
         end
 end
 

%%
for i = 1:6
    itcDiff{i} = itcAllTd{i} - itcAllRett{i}; %differences of averages
end
data = itcDiff; %to plot the difference between the groups

%% run FDR_clustBase_permuTest with the individual structures (itcAllTd, itcAllRett), to get cluster-based permutation test. 
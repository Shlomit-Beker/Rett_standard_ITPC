%% ITPC for Rett participants across all time points. 
% Shlomit Beker shlomitbeker@gmail.com

% Run after Main_Rett.m
% parameters for ITPC
%load 'AllData_long.mat' %need to load the trials that are saved. (and are
%output of Main_Rett. Two options: long (2sec) or short (0.9sec) trials. 
%for collapsed ITPC, better use the short epoch
%AllData = AllData_long; %change if needed
% Shlomit Beker shlomitbeker@gmail.com



if  ~exist('AllData')
    load AllData
else

clear PLallstim PLV groupPLV

LENGTH_WIN = 2;    % length of trial in sec. change according to long (2)/short(0.9)
C = {'FCz';'FC3';'FC4'};
SAMP_RATE = 512; % of the trials, after preprocessing
LOW_FREQUENCY = 2;
HIGH_FREQUENCY = 15; %range of frequencies on which to make the coherence
OMEGA = 6;


TIME_WINDOW = 1:LENGTH_WIN*SAMP_RATE;
CONDNUM = 6;
RAND = 1;
%% number of subjects (ctrl; rett) reduced from 18 to 17 due to bad participant. 
AllData(:,25) = [];
numSub = [24;17]; 

clear i;

labels = cell(0);
for i = 1:size(AllData{1,1}.chanlocs,2)
    labels{i} = AllData{1,1}.chanlocs(i).labels;
end
CHANNELS = [find(strcmp(labels,C{1})),find(strcmp(labels,C{2})),find(strcmp(labels,C{3}))];
cd '/Users/shlomit/Dropbox (EinsteinMed)/TUFI/Rett_ITPC';
cd 'C:\Users\sbeker\Dropbox (EinsteinMed)\TUFI\Rett_ITPC';
%% Create the data mat files, by condition

tic
flag = 0;
participantFlag = 0;
currentDataRand = zeros(0,0,0);
PLV = cell(0);
for COND = 1:CONDNUM % number of different conditions (3ISI*2std/dev)
    STphase_all_stim = [];
        for participant = 1:size(AllData,2)
            clear STphase1
            currentData = AllData{COND,participant}.data; %data from one participant and one condition (isi+dev/st). 
            Ntrial = size(currentData,3);
            A = randperm(Ntrial);
             if RAND
                  for i = 1:size(currentData,3)
                        currentDataRand(:,:,i) = currentData(:,:,A(i)); % to rand order of trial within condition. So correlation wont be affected by the order 
                  end
             end
            
            for i_trial = 1:Ntrial
                    flag = flag+1;
                    [wave,period,scale,cone_of_influence] = basewave4(squeeze(mean(currentData(CHANNELS,:,i_trial)))',SAMP_RATE,LOW_FREQUENCY,HIGH_FREQUENCY,OMEGA,0);
                    STphase1(:,:,i_trial) = squeeze(angle(wave));
            end
            frequencies = 1./period;
            
            PL = squeeze(mean((abs(sum(exp(1i*STphase1(:,TIME_WINDOW,:)),3))/Ntrial),2));
            PLallstim(:) = PL;
            PLV{COND,participant} = PLallstim;
            
        end
end
toc

%% Plot ITPCs for the two groups per condition
clear groupPLV groupPLVerr

% Means and Plot
for j = 1:size(AllData,1)
    meansPLV{j} = cell2mat(PLV(j,:)');
end

%TD
for k = 1:CONDNUM
    subjectsPLV{1}{k} = meansPLV{k}(1:numSub(1),:);
    groupPLV{1}(k,:) = mean(meansPLV{k}(1:numSub(1),:)); %mean
    groupPLVerr{1}(k,:) = std(meansPLV{k}(1:numSub(1),:))./sqrt(numSub(1)); %standard error of the mean
end
 
%RETT
for k = 1:CONDNUM
    subjectsPLV{2}{k} = meansPLV{k}(numSub(1)+1:numSub(2)+numSub(1),:);
    groupPLV{2}(k,:) = mean(meansPLV{k}(numSub(1)+1:numSub(2)+numSub(1),:)); %mean
    groupPLVerr{2}(k,:) = std(meansPLV{k}(numSub(1)+1:numSub(2)+numSub(1),:))./sqrt(numSub(2)); %standard error of the mean
end

COLORS = {'k',[122,122,122]./255};
TITLES = {'ISI 450 std','ISI 450 dev','ISI 900 std','ISI 900 dev','ISI 1800 std','ISI 1800 dev'};

presOrder = [1,3,5,2,4,6];
figure
for k = 1:CONDNUM
    subplot(2,3,k)
    
   % TD group
    
    Am = groupPLV{1}(presOrder(k),:);
    Ae = groupPLVerr{1}(presOrder(k),:);
    [mean_smooth{1}, error_smooth] = drawBoundedLines_NEW(Am,Ae,SAMP_RATE); % calc bounded lines
    shadedErrorBar(frequencies,mean_smooth{1},error_smooth,{COLORS{1},'LineWidth',1},0); %plot them
    
    hold on
    
    % Rett group
    
    Bm = groupPLV{2}(presOrder(k),:);
    Be = groupPLVerr{2}(presOrder(k),:);
    [mean_smooth{2}, error_smooth] = drawBoundedLines_NEW(Bm,Be,SAMP_RATE); % Draw bounded lines
    shadedErrorBar(frequencies,mean_smooth{2},error_smooth,{COLORS{2},'LineWidth',1},1); %plot them
    
    %plotting configurations
    
    hold on;
    text(5.5,0.25,'p < 0.01','Color','r','FontSize',12);
    ylim([0.04 0.3])
    xlim([LOW_FREQUENCY HIGH_FREQUENCY])
    title(TITLES{presOrder(k)})
    xlabel('Frequency');
    ylabel('ITPC');
    set(gca,'fontsize', 14);

end

sgtitle('ITPC per condition')

%% Permutation test - run on each k condition (ISI), and calc stats for permutation test comparison between the groups. 
p = [];
observeddifference = [];
effectsize = [];
for k = 1:CONDNUM
    [p(k), observeddifference(k), effectsize(k)] = permutationTest(groupPLV{1}(presOrder(k),:),groupPLV{2}(presOrder(k),:),1000);
end


%% Plot ITPC for each condition and individual subjects.
TITLES = {'ISI 450 NT','ISI 900 NT','ISI 1800 NT','ISI 450 Rett','ISI 900 Rett','ISI 1800 Rett'};

fig = figure; 
%sgtitle('ITPC - Standard')
suptitle('ITPC - Standard')


flag = 0;
for i = 1:length(subjectsPLV)
    for j = [1,3,5] %std
        flag = flag+1
        subplot(2,3,flag)
        plot(frequencies,subjectsPLV{i}{j})
        title(TITLES{flag})
        xlim([2 12])
        ylim([0 0.5])
        xlabel('Freq (Hz)');
        ylabel('ITPC (AU)');
        hold on;
    end
end




fig = figure; 
%sgtitle('ITPC - Deviant')
suptitle('ITPC - Deviant')

flag = 0;
for i = 1:length(subjectsPLV)
    for j = [2,4,6] %std
        flag = flag+1
        subplot(2,3,flag)
        plot(frequencies,subjectsPLV{i}{j})
        title(TITLES{flag})
        xlim([2 12])
        ylim([0 0.5])
        xlabel('Freq (Hz)');
        ylabel('ITPC (AU)');
        hold on;
    end
end


%% Plot ITPC for all conditions per group. Change GROUP between 1 (TD) and 2 (Rett)

GROUP = 1;

if GROUP == 1
    n = 1:1:numSub(1);
    Title = 'TD group';
else 
    n = numSub(1)+1:numSub(2)+numSub(1);
    Title =  'Rett group';
end

figure;
clear groupPLV
COLORS = {'b','b','r','r','g','g'};
groupPLV = [];
for j = 1:length(meansPLV)
    if mod(j,2) == 0
        width = 2;
        specifier = '--';
    else
        width = 2;
        specifier = '-';
    end
    groupPLV(j,:) = mean(meansPLV{j}(n,:));
    %plot(frequencies,groupPLV(j,:),'Color',COLORS{j},'LineWidth',width);
    plot(frequencies,groupPLV(j,:),specifier, 'Color',COLORS{j},'LineWidth',width);
    hold on; 
end
legend('ISI 450 std','ISI 450 dev','ISI 900 std','ISI 900 dev','ISI 1800 std','ISI 1800 dev');
title(Title);
xlabel('Frequencies');
ylabel('ITPC');
set(gca,'fontsize', 15);
ylim([0 0.3]);
xlim([LOW_FREQUENCY HIGH_FREQUENCY]);




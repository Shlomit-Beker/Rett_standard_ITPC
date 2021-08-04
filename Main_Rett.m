% Main_Rett code to load EEG data from all subjects. creates AllData file. 
% Shlomit Beker shlomitbeker@gmail.com



tic

% Setting path configurations. Change according to paths. 

PATH_data = 'C:\Users\sbeker\Dropbox (EinsteinMed)\TUFI\Rett_ITPC\DATA short'; % Lab computer path - short or long

PATH_data = '/Users/shlomit/Dropbox (EinsteinMed)/TUFI/Rett_ITPC/DATA long';   % Mac path - short or long

addpath('C:\Users\sbeker\Dropbox (EinsteinMed)\GENERAL ANALYSIS FILES');
addpath('/Users/shlomit/Dropbox (EinsteinMed)/GENERAL ANALYSIS FILES');
startup;

%% This part creates the AllData structure and ID list from the individual data structures

cd(PATH_data);
subjectFile = dir(PATH_data);
AllData = [];
subNum = 0;
for k = 4:size(subjectFile,1)-2
    currentSub = subjectFile(k).name;
    load(fullfile(PATH_data,currentSub));
    subNum = subNum + 1;
    AllData = cat(2,AllData,ERPs); % cell array with all subjects in all conditions
    a = regexp(currentSub,'\d');
    Names_long(subNum,1) = str2num(currentSub(a));
end
toc    
     
%% number of trials
for i = 1:size(AllData,1)
    for j = 1:size(AllData,2)
        trialNumber(i,j) = size(AllData{i,j}.data,3);
    end
end

%% in case AllData is not a structure
for i = 1:size(AllData,1)
    for j = 1:size(AllData,2)
        trialNumber(i,j) = size(AllData{i,j},3);
        
    end
end
    
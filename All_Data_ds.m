% downsampling data. 
% input: pre-ds data
% shlomitbeker@gmail.com

function AllData_ds = All_Data_ds(AllData)

for i = 1:size(AllData,1)
    for j = 1:size(AllData,2)
        AllData_ds{i,j} = AllData{i,j}.data(:,1:2:end,:);  %DOWN SAMPLED
    end
end
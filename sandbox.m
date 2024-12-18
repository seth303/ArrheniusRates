% % sandbox

%% Manual Plot
clc; clear; close all;
saveData = 1;
initialGuess = [1e-11 250e3 0.7];
overwriteIC = 0;
saveFMinSearch = 0;
plotIntermediate = 1;
silenceOutput = 0;
plotMultipleEEDF = 0;
reactionToPlot = 132;
load('MASTER_BOLSIG_Rates_Data.mat')

opts = [saveData, reactionToPlot, initialGuess, overwriteIC, saveFMinSearch, plotIntermediate, silenceOutput, plotMultipleEEDF];
[reactionTable, data] = fitReactionsWithDistribution(data, opts);
disp('No distribution fitting')
figure;
[reactionTable, data] = fitReactions(data, opts);

%% Check for NaN's
clc; clear; close all;
saveData = 1;
initialGuess = [1e-11 250e3 0.7];
overwriteIC = 1;
saveFMinSearch = 0;
plotIntermediate = 0;
silenceOutput = 0;
plotMultipleEEDF = 0;
reactionToPlot = 209;
load('MASTER_BOLSIG_Rates_Data.mat')

opts = [saveData, reactionToPlot, initialGuess, overwriteIC, saveFMinSearch, plotIntermediate, silenceOutput, plotMultipleEEDF];
[reactionTable, data] = fitReactions(data, opts);

save('MASTER_reactionTable.mat',"reactionTable")
save('MASTER_BOLSIG_Rates_Data.mat',"data")

nanIdX = isnan(reactionTable.A);
nanTable = reactionTable(nanIdX, :)

%% Check for Duplicates
load('MASTER_BOLSIG_Rates_Data.mat')
clc
chemistryData = cell(length(fieldnames(data)),1);
for i = 1:length(fieldnames(data))
    reactionField = sprintf('Reaction_%d', i);
    chemistryData{i} = data.(reactionField).chemistry;
end

% check for duplicates index
[U,I] = unique(chemistryData, 'first');
x = 1:length(chemistryData);
x(I) = []
load MASTER_reactionTable.mat
n = 1;
searchString = table2array(reactionTable(x(n),2));

reactionTable(matches(reactionTable.Chemistry, searchString),:)

%% remove duplicates in the struct
% Assume `reactions` is your structure
% Example duplicate reactions to remove
duplicateIndices = 34; % Indices of duplicates to remove

% Remove duplicates
data = rmfield(data, strcat("Reaction_", string(duplicateIndices)));

% Get the remaining field names and sort by numerical index
remainingFields = fieldnames(data);
numFields = length(remainingFields);

% Initialize a new structure to hold reordered fields
newReactions = struct();

% Loop to rename fields sequentially
for i = 1:numFields
    % Create new field name
    newFieldName = strcat("Reaction_", string(i));
    % Copy data to new structure with sequential field names
    newReactions.(newFieldName) = data.(remainingFields{i});
end

% Replace original structure with reordered structure
data = newReactions;
save('MASTER_BOLSIG_Rates_Data.mat',"data")

%% check for shortfields
fieldLength = zeros(length(fieldnames(combinedData)),1);
for i = 1:length(fieldnames(combinedData))
    reactionField = sprintf('Reaction_%d', i);
    fieldLength(i) = size(fieldnames(combinedData.(reactionField)),1);
end

shortFields = find(fieldLength == 3)


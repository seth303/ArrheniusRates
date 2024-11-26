function [matchingChemistry, idx1, idx2, updatedStruct2, missingIdx] = compareChemistryFields(struct2, forcedFit)
% load the database
load 'MASTER_BOLSIG_Rates_Data.mat' data
struct1 = data;

% Get the number of reactions in both structures
numReactions1 = numel(fieldnames(struct1));
numReactions2 = numel(fieldnames(struct2));

% Initialize arrays to store Chemistry values
chemArray1 = cell(numReactions1, 1);
chemArray2 = cell(numReactions2, 1);

% Extract Chemistry values from struct1
for i = 1:numReactions1
    reactionField1 = sprintf('Reaction_%d', i);
    chemArray1{i} = struct1.(reactionField1).chemistry;
end


% Extract Chemistry values from struct2
for i = 1:numReactions2
    reactionField2 = sprintf('Reaction_%d', i);
    chemArray2{i} = struct2.(reactionField2).chemistry;
end
% Find matching Chemistry strings
[matchingChemistry, idx1, idx2] = intersect(chemArray1, chemArray2);
disp(['Found ', num2str(length(matchingChemistry)), ' matching reactions out of ' , num2str(length(fields(struct2)))])

if length(matchingChemistry) ~= length(fields(struct2))
    missingIdx = find(~matches(chemArray2,chemArray1));
else
    missingIdx = 0;
end


% Initialize updatedStruct2 as a copy of struct2
updatedStruct2 = struct2;

% Copy data from struct1 to struct2 for matching Chemistry fields
for i = 1:numel(idx1)
    reactionField1 = sprintf('Reaction_%d', idx1(i));
    reactionField2 = sprintf('Reaction_%d', idx2(i));

    % Copy data from struct1 to the corresponding field in struct2
    updatedStruct2.(reactionField2) = struct1.(reactionField1);
end
% also want to copy over noncombined data into the struct
copiedFields = 0;
for i = 1:numel(chemArray1)
    % for all the values that are not in idx1, that is the values that are
    % only present in struct1, copy them over at the end of updateStruct2

    % ex, chemArray1(1) = N2O -> N2O+ which is the same as chemArray2(13)
    % so this needs to be skipped because we already copied it over above
    % but chemArray1(2) = NO2 which is not in chemArray2 so that one needs
    % to be copied to the end of updateStruct2
    if ~any(idx1==i)
        copiedFields = copiedFields + 1;

        % if any value of idx1 is NOT equal to the current iteration, copy
        % that data from struct1 to the end of struct2
        reactionField1 = sprintf('Reaction_%d', i);
        reactionField2 = sprintf('Reaction_%d', length(chemArray2)+copiedFields);
        updatedStruct2.(reactionField2) = struct1.(reactionField1);
    end
end
if forcedFit
    % if forcedChemistry flag is true, set all fields in chemArray2 to 'none'
    [chemArray1{:}] = deal('chemArray1');
    [chemArray2{:}] = deal('chemArray2');
    missingIdx = 1:length(fieldnames(struct1));
    disp(['Forced Fit Initiated. Fitting ', num2str(length(fieldnames(updatedStruct2))), ' reactions'])
end
end
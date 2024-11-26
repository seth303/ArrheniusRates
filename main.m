%% Step 0: Restart
clc; clear; close all;
%% Step 1: import cross sectional data from input bolsig file
newData = importCrossSectionData('bolsigdb.dat');
%% Step 2: check if any of the imported data matches already fitted data
% if there are matches, import the data into newData to form combinedData
% if I want to force the fitting algorithm or all data
forcedFit = 0;
[matchingChemistry, idx1, idx2, combinedData, missingIdx] = compareChemistryFields(newData, forcedFit);
%% Step 3: Check If New Data Exists
newReactionsToFit = checkForNewReactions(matchingChemistry,newData,combinedData,missingIdx);
%% Step 3.1: Choose which reactions to review (optional)
reviewFit = 0;
startIdx = 19;
endIdx = 20;
load("MASTER_reactionTable.mat");
[C, idxA,idxB,bolsigRxnsName] =  findBolsigRxns(reactionTable);
% sort
[idxB,sortIdx] = sort(idxB,'ascend');
idxA = idxA(sortIdx);
C = C(sortIdx);
%% Step 4: If so, Fit it
saveData = 0;
if newReactionsToFit || forcedFit || reviewFit
    % options
    initialGuess = [9.71e-15 17e3 1.02];
    overwriteIC = 0;
    saveFMinSearch = 0;
    plotIntermediate = 0;
    silenceOutput = 0;
    plotMultipleEEDF = 0;
    fitType = 'inversegaussian';
    plotTrapz = 0;
    shiftFactor = 0;
    load('validFitTypes.mat')
    upperLimit = 100;

    if ~ismember(fitType, validFitTypes) 
        disp('Valid Fit Types: ')
        disp(string(validFitTypes')); 
        error('Invalid Fit Type.'); 
    else
        validFitTypeChosen = 1;
    end

    if forcedFit && saveData
        warning('Fitting algorithm forced and data will be saved. Ensure all data is backed up.')
        input('Would you like to continue?')
    end

    if reviewFit
        disp('Review Fit: User has defined reactions to fit: ')
        disp(C)
        missingIdx = idxA;
    end

    % for i = startIdx:length(missingIdx)
    for i = startIdx:endIdx
        close all;
        % assign the missing reaction to "reactionToPlot"
        reactionToPlot = missingIdx(i);
        % define opts
        opts = {saveData, reactionToPlot, initialGuess, overwriteIC, saveFMinSearch, plotIntermediate, silenceOutput, plotMultipleEEDF,fitType,plotTrapz,shiftFactor,upperLimit};
        [reactionTable, combinedData, exitCondition] = fitReactionsWithSmoothingSpline(combinedData, opts);
       
        % check for user approval
        [combinedData,reactionTable] = userApproval(combinedData,i,opts);
       
    end
else
    disp('No new reactions to fit, moving on to kinetics file editing')
end
%% Step 5, Update the Database
if newReactionsToFit || saveData
    updateDataBase(combinedData,reactionTable, missingIdx)
else
    disp('No new reactions to fit, kinetics file generated')
end

%% Step 6, Modify the Kinetics File
originalFile = 'kin1NxO.inp';
copiedFile = 'mod_kin1NxO.inp';

refactoringKineticsFile(originalFile, copiedFile)
disp('End Program')
% close all
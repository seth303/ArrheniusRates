function [combinedData,reactionTable] = userApproval(fittedData,i,opts)
%% Options Assignment
saveData = opts{1};
reactionToPlot = opts{2};
initialGuess = opts{3};
overwriteIC = opts{4};
saveFMinSearch = opts{5};

plotIntermediate = opts{6};
silenceOutput = opts{7};
plotMultipleEEDF = opts{8};
fitType = opts{9};
plotTrapz = opts{10};
shiftFactor = opts{11};
upperLimit = opts{12};
load 'validFitTypes.mat' validFitTypes
%%
proceed = input('Press Enter to accept fit, or any other key to provide new initial conditions: ', 's');
if isempty(proceed)  % If Enter is pressed (empty input)
    disp(['Continuing the loop, iteration: ', num2str(i)]);
    % need to add the recently created fit to the combinedData
    combinedData = fittedData;
    reactionTable = manualReactionTable(combinedData);
else
    userApproved = false;
    % plotIntermediate = 1;
    % opts(6) = {1};
    while(~userApproved)
        newInitial_A = input(['Provide new value of A. Previous Value: ',num2str(initialGuess(1)),'. New Value: ']);
        newInitial_Ea = input(['Provide new value of Ea. Previous Value: ',num2str(initialGuess(2)),'. New Value: ']);
        newInitial_n = input(['Provide new value of n. Previous Value: ',num2str(initialGuess(3)),'. New Value: ']);
        % newFitType = input(['Provide new fit type. Previous Value: ',fitType,'. New Value: '],"s");
        % newShiftFactor = input(['Provide new shift factor. Previous Value: ',num2str(shiftFactor),'. New Value: ']);
        newLowerLimit = input(['Provide new lower limit. Previous Value: ',num2str(upperLimit),'. New Value: ']);
        if ~isempty(newInitial_A); initialGuess(1) = newInitial_A; end
        if ~isempty(newInitial_Ea); initialGuess(2) = newInitial_Ea; end
        if ~isempty(newInitial_n); initialGuess(3) = newInitial_n; end
        % if ~isempty(newFitType); fitType = newFitType; end
        % if ~isempty(newShiftFactor); shiftFactor = newShiftFactor; end
        if ~isempty(newLowerLimit); upperLimit = newLowerLimit; end

        if ~ismember(fitType, validFitTypes)
            disp('Valid Fit Types: ')
            disp(string(validFitTypes'));
            validFitTypeChosen = 0;
        else
            validFitTypeChosen = 1;
        end

        while (~validFitTypeChosen)
            disp('Invalid Fit Type Chosen')
            newFitType = input(['Provide new fit type. Previous Value: ',fitType,'. New Value: '],"s");
            if ~isempty(newFitType); fitType = newFitType; end
            if ~ismember(fitType, validFitTypes)
                disp('Valid Fit Types: ')
                disp(string(validFitTypes'));
                validFitTypeChosen = 0;
            else
                validFitTypeChosen = 1;
            end
        end

        close all;
        overwriteIC = 1;
        opts = {saveData, reactionToPlot, initialGuess, overwriteIC, saveFMinSearch, plotIntermediate, silenceOutput, plotMultipleEEDF,fitType,plotTrapz,shiftFactor,upperLimit};
        [reactionTable, fittedData, exitCondition] = fitReactionsWithSmoothingSpline(fittedData, opts);
        approval = input('Fit Accepted? [y/n]: ', 's');

        if strcmp(approval, 'y')
            % plotIntermediate = 0;
            userApproved = true;
            % overwriteIC = 0;
            initialGuess = [1e-11 180e3 0.3];
            % need to add the recently created fit to the combinedData
            combinedData = fittedData;
            break;
        else
            userApproved = false;
        end
    end
end
end
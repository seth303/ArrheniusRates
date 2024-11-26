function [newReactionsToFit] = checkForNewReactions(matchingChemistry,newData,combinedData,missingIdx)

if length(matchingChemistry) == length(fields(newData))
    newReactionsToFit = false;
else
    newReactionsToFit = true;
    disp('Missing Fits for: ');
    for i = 1:length(missingIdx)
        reactionField_Missing = sprintf('Reaction_%d', missingIdx(i));
        disp(combinedData.(reactionField_Missing).chemistry);
    end
end

end
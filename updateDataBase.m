function updateDataBase(combinedData,reactionTable, missingIdx)

disp('Updating Database with reactions for: ')
for i = 1:length(missingIdx)
    if missingIdx ~= 0
        reactionField_Missing = sprintf('Reaction_%d', missingIdx(i));
        disp(combinedData.(reactionField_Missing).chemistry);
    else
        disp('No New Reactions')
    end
end
save('MASTER_reactionTable.mat',"reactionTable")
data = combinedData;
save('MASTER_BOLSIG_Rates_Data.mat',"data")

end
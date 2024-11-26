function refactoringKineticsFile(originalFile, copiedFile)
load MASTER_reactionTable.mat reactionTable
removeBolsigReactionGroups(originalFile, copiedFile, 0);
reactionsReplaced =  1;
maxReplacement = 999;
%% Find matching reactions
fprintf(' \n')
fileID = fopen(copiedFile, 'r');  % Open the copied file for reading
tempFileID = fopen('tempfile.txt', 'w');  % Temporary file for writing modified content
currLinNum = 0;
searchStringArray = (reactionTable.Chemistry);

% for i = 1:size(reactionTable.Chemistry,1)
%     % Define the search string you are looking for
%     searchString = char(reactionTable.Chemistry(i));
%     % disp(['Searching for Reaction_', num2str(i), ': ', searchString])
%     % Step 3: Process the copied file line by line
while ~feof(fileID)
    % Read the current line
    currentLine = fgetl(fileID);
    currLinNum = currLinNum + 1;

    % I want to check if the current line contains any of the
    % reactions in the table, and if so, pick the first one and use it
    % for the replacement
    [isMatch, idx] = searchLineForPartialMatch(currentLine, searchStringArray);

    if (maxReplacement < reactionsReplaced)
        disp('Maximum reactions replaced')
        isMatch = false;
    end

    if isMatch && ~contains(currentLine, '#') && contains(currentLine, 'BOLSIG')
        searchString = char(reactionTable.Chemistry(idx));
        matchOutput = ['Found matching Reaction_', num2str(idx), ': ', searchString, ' on line ',  num2str(currLinNum+1) ];
        disp(matchOutput)

        fprintf(' \n')
        % fprintf(tempFileID, '# %s\n', matchOutput);
        % now that we have found a match, we need to replace the '!
        % BOLSIG XXXX' into the arrenhius form of the equation using
        % the values from the reactionTable
        %  ex: if A = 11.853, Ea = 1.84E5 and n = -3.27 then we would
        %  want to the new reaction rate in the kinetics file to be in
        %  the form of:
        % ' ! 11.853D0*Te^(-3.27D0)*exp(-1.84D5/Te)
        % Split the current line at '!' (before ! is the base part)
        parts = strsplit(currentLine, '!');
        baseLine = (parts{1});
        % custom function that converts numbers to strings in ZDPlasKin
        % compatible format
        str_A = toScientificNotation(reactionTable.A(idx));
        str_Ea = toScientificNotation(reactionTable.Ea(idx));
        str_n = toScientificNotation(reactionTable.n(idx));
        chemistry = char(reactionTable.Chemistry(idx));
        % construct the newline with reaction rates
        % in Fortran, ** is the same as ^
        % ^ will not work
        newLine = [baseLine, '! ', str_A, '*Te**(', str_n,')*exp(-',str_Ea,'/Te)'];
        newLineChemistry = ['# Previous Reaction: ',chemistry];
        % Write each expanded line to the temp file
        fprintf(tempFileID, '%s\n', newLineChemistry);
        fprintf(tempFileID, '%s\n', newLine);
        currLinNum = currLinNum + 1;
        % if this is true, we need to rewind the line loader, and
        % % iterate to the next reaction in the table
        % frewind(fileID)
        % currLinNum = 0;
        % break
        reactionsReplaced =  reactionsReplaced +1;
    else
        % If the search string is not found, write the current line as is to the temp file
        fprintf(tempFileID, '%s\n', currentLine);
    end
end

% % if no match is found, rewind the line loader
% frewind(fileID)
% currLinNum = 0;

% end
% Step 4: Close both files
fclose(fileID);
fclose(tempFileID);

% Step 5: Replace the original copy with the modified temporary file
movefile('tempfile.txt', copiedFile);  % Overwrite the copied file with the modified content

disp(['Modifications complete! The modified file is saved as ', copiedFile]);
end
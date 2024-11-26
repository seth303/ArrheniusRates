function removeBolsigReactionGroups(originalFile,copiedFile, outputFlag)
% Define the original and copied file names


currLinNum = 0;

% Step 1: Create a copy of the original file
copyfile(originalFile, copiedFile);

% Step 2: Open the copied file for reading and also open a temporary file for writing modifications
fileID = fopen(copiedFile, 'r');  % Open the copied file for reading
tempFileID = fopen('tempfile.txt', 'w');  % Temporary file for writing modified content

% Define the search string you are looking for
searchString = 'BOLSIG';

% Step 3: Process the copied file line by line
while ~feof(fileID)
    % Read the current line
    currLinNum = currLinNum + 1;
    currentLine = fgetl(fileID);


    % Check if the current line contains the search string
    if contains(currentLine, searchString) && ~contains(currentLine, '#') && contains(currentLine, '!')
        % Placeholder for your action when the string is found
        if outputFlag; disp(['String found on line ', num2str(currLinNum), ': ', currentLine]); end
        % Further check if the line contains '@R'
        if contains(currentLine, '@R')
            if contains(currentLine, '@B')
                if outputFlag; disp(['Double group found on line ', num2str(currLinNum), ': ', currentLine]); end
                % if a B is found, we want to modify the baseLine to
                % include the B reactions similar to how we are modifying
                % the R line.
                fprintf(tempFileID, '# The following section was modified to ungroup the reactions\n');
                parts = strsplit(currentLine, '@B');
                baseLine = strtrim(parts{1});
                reactionLine = strsplit(parts{2}, '@R');

                % Read the next line, which contains the variations
                nextLine = fgetl(fileID);
                variations_B = strtrim(nextLine); % Variations are on the next line

                % Read the line after that to get the @R values
                nextLine = fgetl(fileID);
                variations_R = strtrim(nextLine); % Variations are on the next line

                % Split the variations into separate entries
                variationList_B = strsplit(variations_B);
                variationList_R = strsplit(variations_R);

                % Write the expanded lines to the temporary file
                for i = 3:length(variationList_B)
                    newLine = [baseLine, ' ', variationList_B{i}, reactionLine{1}, variationList_R{i}, ]; % Combine the base line with each variation
                    fprintf(tempFileID, '%s\n', newLine);        % Write each expanded line to the temp file

                end
            else % meaning we found a @R but not @B

                fprintf(tempFileID, '# The following section was modified to ungroup the reactions\n');
                % Split the current line at '@R' (before @R is the base part)
                parts = strsplit(currentLine, '@R');
                baseLine = strtrim(parts{1});

                % Read the next line, which contains the variations
                nextLine = fgetl(fileID);
                variations = strtrim(nextLine); % Variations are on the next line

                % Split the variations into separate entries
                variationList = strsplit(variations);

                % Write the expanded lines to the temporary file
                for i = 3:length(variationList)
                    newLine = [baseLine, ' ', variationList{i}]; % Combine the base line with each variation
                    fprintf(tempFileID, '%s\n', newLine);        % Write each expanded line to the temp file

                end
            end
        else
            % If no '@R' is found, just write the current line as is to the temp file
            fprintf(tempFileID, '%s\n', currentLine);
        end
    else
        % If the search string is not found, write the current line as is to the temp file
        fprintf(tempFileID, '%s\n', currentLine);
    end
end

% Step 4: Close both files
fclose(fileID);
fclose(tempFileID);

% Step 5: Replace the original copy with the modified temporary file
movefile('tempfile.txt', copiedFile);  % Overwrite the copied file with the modified content

% disp(['Modifications complete! The modified file is saved as ', copiedFile]);
end
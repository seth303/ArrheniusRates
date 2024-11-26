function reactions = importCrossSectionData(filename)
    % Open the text file
    fid = fopen(filename, 'r');
    
    reactions = struct(); % Initialize an empty structure
    reactionCounter = 1;  % Counter for naming reactions (Reaction_1, Reaction_2, etc.)
    
    while ~feof(fid)
        % Read lines and look for "EXCITATION", "ATTACHMENT", or "IONIZATION"
        line = fgetl(fid);
        
        if (contains(line, 'EXCITATION') || contains(line, 'ATTACHMENT') || contains(line, 'IONIZATION')) && (~contains(line, '#'))
            % Store the method (excitation, attachment, or ionization)
            if contains(line, 'EXCITATION')
                method = 'EXCITATION';
            elseif contains(line, 'ATTACHMENT')
                method = 'ATTACHMENT';
            elseif contains(line, 'IONIZATION')
                method = 'IONIZATION';
            end
            
            % The reaction name is the line immediately after the method line
            reactionNameLine = fgetl(fid);
            reactionName = strtrim(reactionNameLine); % Keep the reaction as is, without any modifications
            
            % Skip until we reach the dashed line (beginning of data)
            while ~feof(fid)
                line = fgetl(fid);
                if contains(line, '-----------------------------')
                    break;
                end
            end
            
            % Read the numeric data (after the dashed line)
            data = textscan(fid, '%f %f', 'Delimiter', '\t', 'CollectOutput', true);
            
            % Store the numeric matrix and the original chemistry string
            numericMatrix = data{1}; % Extract the numeric matrix from the cell
            reactionFieldName = sprintf('Reaction_%d', reactionCounter); % Create Reaction_1, Reaction_2, etc.
            
            % Assign the data, chemistry, and method to the structure
            reactions.(reactionFieldName).data = numericMatrix;
            reactions.(reactionFieldName).chemistry = reactionName; % Store the unmodified reaction string
            reactions.(reactionFieldName).method = method;   % Store the method (excitation, attachment, ionization)
            
            % Increment the counter for the next reaction
            reactionCounter = reactionCounter + 1;
        end
    end
    
    % Close the file
    fclose(fid);
end

function duplicates = findDuplicateSubstrings(inputStr, minLength)
    % Initialize output as an empty cell array
    duplicates = {};
    n = strlength(inputStr);
    
    % Loop over possible substring lengths (from minLength to half of the string length)
    for len = minLength:floor(n/2)
        for i = 1:(n - len + 1)
            % Extract the current substring
            substring = extractBetween(inputStr, i, i + len - 1);
            
            % Find all occurrences of this substring in the string
            matches = regexp(inputStr, substring, 'start');
            
            % If there are more than one occurrence, store the substring and positions
            if numel(matches) > 1
                duplicates{end+1, 1} = substring;   % Store the substring
                duplicates{end, 2} = matches;       % Store starting positions of each occurrence
            end
        end
    end
end
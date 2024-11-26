function updatedArray = addSuffixToDuplicates(cellArray)
    % Initialize output array
    updatedArray = cellArray;
    
    % Create a map to track counts of each name in column 2
    countMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    for i = 1:size(cellArray, 1)
        entry = cellArray{i, 2};
        
        % Check if the entry already exists in the map
        if isKey(countMap, entry)
            % Increment the count and add it as a suffix
            countMap(entry) = countMap(entry) + 1;
            suffix = countMap(entry);
            updatedArray{i, 2} = sprintf('%s_%d', entry, suffix);
        else
            % Add the entry to the map with a count of 1
            countMap(entry) = 1;
            updatedArray{i, 2} = entry;
        end
    end
end
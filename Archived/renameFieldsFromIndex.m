function structOut = renameFieldsFromIndex(structIn, startIdx)
    % Get all the field names
    fieldNames = fieldnames(structIn);
    
    % Loop through fields starting from 'startIdx' down to the last field
    for i = numel(fieldNames)+1:-1:startIdx
        % Generate old and new field names
        oldFieldName = sprintf('Reaction_%d', i);
        newFieldName = sprintf('Reaction_%d', i - 1);
        
        % Check if the old field exists
        if isfield(structIn, oldFieldName)
            % Rename the field: copy the data to the new field and remove the old one
            structIn.(newFieldName) = structIn.(oldFieldName);
            structIn = rmfield(structIn, oldFieldName);
        end
    end
    
    % Output the updated structure
    structOut = structIn;
end

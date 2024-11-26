function [reactionTable] = manualReactionTable(data)
%% Create Reaction Table
% Assume 'data' is the input 1x1 structure with 67 fields.
% Each field has subfields '.chemistry', '.A', '.Ea', and '.n'.

% Get all field names of the structure
fields = fieldnames(data);

% Preallocate cell arrays for the table columns
fieldNames = cell(length(fields), 1);
chemistryValues = cell(length(fields), 1);
A_values = zeros(length(fields), 1);  % Assuming numeric values for A
Ea_values = zeros(length(fields), 1); % Assuming numeric values for Ea
n_values = zeros(length(fields), 1);  % Assuming numeric values for n

% Loop over all fields and extract the relevant values
for i = 1:length(fields)
    fieldName = fields{i};
    fieldNames{i} = fieldName;  % Store the field name
    % if the fitting exists, that is, if the data.(fieldName) length is
    % greater than 3, save the data, else, input NaN
    if length(fieldnames(data.(fieldName))) > 3
        chemistryValues{i} = data.(fieldName).chemistry;  % Store the chemistry value
        A_values(i) = data.(fieldName).A;  % Store the A value
        Ea_values(i) = data.(fieldName).Ea;  % Store the Ea value
        n_values(i) = data.(fieldName).n;  % Store the n value
    else
        chemistryValues{i} = data.(fieldName).chemistry;  % Store the chemistry value
        A_values(i) = NaN;  % Store the A value as NaN if no fit exists
        Ea_values(i) = NaN;  % Store the Ea value as NaN if no fit exists
        n_values(i) = NaN;  % Store the n value as NaN if no fit exists
    end
end

% Create a table with five columns: 'ReactionName', 'Chemistry', 'A', 'Ea', and 'n'
reactionTable = table(fieldNames, chemistryValues, A_values, Ea_values, n_values, ...
    'VariableNames', {'ReactionName', 'Chemistry', 'A', 'Ea', 'n'});

% Display the table
% disp(reactionTable);
save("MASTER_reactionTable",'reactionTable')
end
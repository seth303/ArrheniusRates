function [isMatch, idx] = searchLineForPartialMatch(line, cellArray)
    % Initialize output variables
    % isMatch = false;
    % idx = 0;
    line = strsplit(line, 'BOLSIG');
    if length(line) > 1
        line = strtrim(line{2});
        % Use regular expression to replace multiple spaces with a single space
        line = regexprep(line, '\s+', ' ');
    end

    % alternative way to look for matches
    idx = find(ismember(string(cellArray),string(line)));
    if ~isempty(idx)
        isMatch = true;
    else
        idx = 0;
        isMatch = false;
    end

    % % Loop through each element in the cellArray to check for partial matches
    % for i = 1:length(cellArray)
    %     % Use regular expression to replace multiple spaces with a single space
    %     cellArray{i} = regexprep(cellArray{i}, '\s+', ' ');
    %     if contains(line, cellArray{i})
    %         isMatch = true;
    %         idx = i;
    %         break;  % Exit the loop once a match is found
    %     end
    % end

end

clc
for i = 23:67
    reactionField = strcat(['Reaction_',num2str(i)]);
    if ~isfield(data.(reactionField), 'SSE')
        i
        error('Data not saved')
    end
end
function [C, idxA,idxB,bolsigRxnsName] =  findBolsigRxns(reactionTable)
% define filepaths
reactionListFilepath = 'Raw QTPlasKin Analysis/ZDPlasKin_FS_Simulations/qt_reactions_list_bolsig.txt';

% import reaction names
varNames = readcell(reactionListFilepath);
varNames = [{0,'Time [s]'};varNames];


% find bolsig reactions
bolsigIdx = contains(string((varNames(:,2))), 'bolsig');
rxnNames = string((varNames));
bolsigRxnsNum = rxnNames(bolsigIdx);
bolsigRxnsName = rxnNames(str2double(bolsigRxnsNum)+1,2);

% split by bolsig
bolsigRxnsName = split(bolsigRxnsName, "bolsig:");
bolsigRxnsName = bolsigRxnsName(:,2);

% define search string
searchArray  = upper(regexprep(reactionTable.Chemistry, '\s+', ''));

%find intersect
[C, idxA, idxB] = intersect(searchArray, bolsigRxnsName);
end
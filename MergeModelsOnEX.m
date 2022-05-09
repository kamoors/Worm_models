function [MergedModel] = MergeModelsOnEX(TargetModel,AddedModel)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


MergedModel = TargetModel;

for n = 1:length(AddedModel.rxns)
    
EX_added_lb = AddedModel.lb(findRxnIDs(AddedModel,AddedModel.rxns(n)));  % Find rxn lb
EX_added_ub = AddedModel.ub(findRxnIDs(AddedModel,AddedModel.rxns(n)));  % Find rxn ub


EX_added_subSystems = AddedModel.subSystems(findRxnIDs(AddedModel,AddedModel.rxns(n))); % Find rxn subSystem
EX_added_rxnName = AddedModel.rxnNames(findRxnIDs(AddedModel,AddedModel.rxns(n))); % Find rxnName
    [x,y] = findMetsFromRxns(AddedModel,char(AddedModel.rxns(n)));
    MergedModel = addReaction(MergedModel,char(AddedModel.rxns(n)), ...
        'metaboliteList',cellstr(x{1}), ...
    'stoichCoeffList',y{1}, ...
    'lowerBound', EX_added_lb, ...
    'upperBound', EX_added_ub, ...
    "subSystem", EX_added_subSystems{1}, ...
    "reactionName", EX_added_rxnName{1});
    %); 
    %, 

    
end



end
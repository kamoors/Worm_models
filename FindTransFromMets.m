function [transRxns, transRxnsBool] = FindTransFromMets(model, metList, targetComps, inclExc, rxnInds, inclObjAsExc, irrevFlag)
% FindTransFromMets 
% Function to find the transport reactions (from one compartment to
% another) from a list of metabolites. 
%   INPUT:
%    model:            COBRA model structure
%
% OPTIONAL INPUT:
%    inclExc:          includes exchange reactions as transport = true
%                      (Default = false)
%    rxnInds:          indices of reactions to test for transport activity.
%                      (default = test all columns of `model.S`)
%    inclObjAsExc:     includes objective as an exchange reaction = true - this is
%                      passed to `findExcRxns`. (default = false)
%    irrevFlag:        model is irreversible = true - this is passed to `findExcRxns`.
%                      (default=false)
%
% OUTPUTS:
%    transRxns:        transport reactions in `rxnInds`
%    nonTransRxns:     non-transport reactions in `rxnInds`
%    transRxnsBool:    checks if `inclExc` or `isNonexchangeTransport` rxns vector is equal to 1
%
% .. Authors:
%       - Jeff Orth, 8/31/07 right now, this function only works with models the compartments [c], [p], and [e].
%       - Jonathan Dreyfuss, 10/9/12modified the function to work with arbitrary compartments, to accept  rxnInds, & to use findExcRxns.
%       - Thierry Mondeel, 07/15/15 modified to also output a boolean vector version of transRxns.
%       - Karlis Moors, 2022 added arguments section; added ability to search from metabolites and specific compartment transport
%       
arguments

model;
metList = model.mets;
targetComps = array2table(model.comps);
inclExc = false;
rxnInds = 1:size(model.S, 2);
inclObjAsExc = false;
irrevFlag = false;
end
















%% Use original function if metList not provided
if isequal(metList, model.mets)
if inclExc
    % findExcRxns returns boolean vector
    isExc0 = findExcRxns(model,inclObjAsExc,irrevFlag);
    % subset to rxnInds
    isExc = isExc0(rxnInds);
else
    isExc=zeros(length(rxnInds), 1);
end

% initialize isNonexchangeTransport rxns vector
isNonexchTrans = zeros(length(rxnInds), 1);
% get compartment symbols for each metabolite
[baseMetNames,compSymbols]=arrayfun(@parseMetNames, model.mets);
for i = 1:length(rxnInds)
    rxnIndTmp=rxnInds(i);
    % get compartment symbols for each rxn
    compSymbolsTmp=compSymbols(model.S(:,rxnIndTmp)~=0);
    % if there's more than 1 compartment involved, it's a transport rxn
    if length(unique(compSymbolsTmp))>1
        isNonexchTrans(i) = 1;
    end
end

% get rxn abbreviations for all rxns in rxnInds
rxnAbbrevs=model.rxns(rxnInds);
% if inclExc==1, exchange rxns will have isExc==1, and should be counted as
% transport rxns; else, all isExc will be 0.
transRxnsBool = isNonexchTrans==1 | isExc==1;
transRxns = rxnAbbrevs(transRxnsBool);
nonTransRxns = setdiff(rxnAbbrevs, transRxns);

return
end
%% If function wants to find from metabolites
if ~isequal(metList, model.mets)

rxnInds = cell([length(metList) 1]);

targetComps = targetComps';
output = struct();

% initialize isNonexchangeTransport rxns vector
isNonexchTrans = cell([length(metList) 1]);

rxnAbbrevs = cell([length(metList) 2]);


transRxns = cell([length(metList) 2]);
%hastargetComps = cell([length(metList) 1]);

for i = 1:length(metList)


[rxns_from_mets, rxnFormulaList] = findRxnsFromMets(model, metList(i));
rxnInds{i} = findRxnIDs(model,rxns_from_mets);



    if inclExc
    % findExcRxns returns boolean vector
    isExc0 = findExcRxns(model,inclObjAsExc,irrevFlag);
    % subset to rxnInds
    isExc = isExc0(rxnInds{i});
    else
    isExc=zeros(length(rxnInds{i}), 1);
    end








% get compartment symbols for each metabolite
[baseMetNames,compSymbols]=arrayfun(@parseMetNames, model.mets);




            for l = 1:length(rxnInds{i})

                    rxnIndTmp=rxnInds{i}(l);

                    % get compartment symbols for each rxn
                    compSymbolsTmp=compSymbols(model.S(:,rxnIndTmp)~=0);
                    % if there's more than 1 compartment involved, it's a transport rxn
                                    if length(unique(compSymbolsTmp))>1
                                        isNonexchTrans{i}(l) = 1;
                                    else
                                        isNonexchTrans{i}(l) = 0;
                                    end
            

                    x = unique(compSymbolsTmp);
                    target_sum = 0;

                        for j = 1:length(targetComps)


                                    for k = 1:length(x)
                                    y = ismember(targetComps{j},x{k});
                                
                                    target_sum = target_sum+y;
                                
                                 
                                    end

                        end


                                    if target_sum == length(targetComps)

                                        isNonexchTrans{i}(l) = 1;
                                
                                    else 

                                        isNonexchTrans{i}(l) = 0;
                                    end

    
            end
rxnAbbrevs{i}=model.rxns(rxnInds{i});
transRxnsBool = isNonexchTrans{i}==1; %| isExc==1;
transRxns{i,1} = rxnAbbrevs{i}(transRxnsBool);
transRxns{i,2} = metList(i);
%
end

for i = 1:length(transRxns)

    if ~isempty(transRxns{i,1})
transRxns{i,1}(:,2) = printRxnFormula(model, transRxns{i});
    end
end


end
end
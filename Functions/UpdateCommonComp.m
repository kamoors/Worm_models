function [updatedModel, dataTable] = UpdateCommonComp(model,OG_mets, Rxns, newCompartmentName, add_COMMON_EX, lb, ub, COMMON_mets)
%UpdateCommonComp
% Function to update the exchange metabolite reactions to create
% COMMON metabolites or update exchange reactions to transport external
% metabolites into a common compartment
%

% model - COBRA model
% OG_mets - original metabolite names met[comp] to be added to [i] (common
% compartment
% newCompartmentName

arguments
    model;
    OG_mets;
    Rxns;
    newCompartmentName = '[i]';
    add_COMMON_EX  = 0  ;
    lb = -1000;
    ub = 1000;
    COMMON_mets = array2table(append(erase(table2array(array2table(erase(table2array(OG_mets), "[e]"))), "[c]"), newCompartmentName));
    
end


updatedModel = model;




%% If no COMMON_mets are added, execute code

if strcmp(COMMON_mets, array2table(append(erase(table2array(array2table(erase(table2array(OG_mets), "[e]"))), "[c]"), newCompartmentName)))   
    dataTable = OG_mets;  % Add iCEL met names


    % for i = 1:height(dataTable) % remove compartment
    % x = strsplit(dataTable(i,1),{'[e0]', '[c0]', '[c]', '[e]', '[m]'});
    % dataTable(i,2) = x(1);
    % end


    dataTable(:,2) = array2table(erase(table2array(dataTable(:,1)), "[e]")); % remove compartment
    dataTable(:,2) = array2table(erase(table2array(dataTable(:,2)), "[c]")); % remove compartment
    dataTable(:,3) = array2table(append(table2array(dataTable(:,2)), newCompartmentName)); %add COMMON compartment
    dataTable(:,4) = Rxns; % get iCEL1314 reaction IDs

    if add_COMMON_EX == 0
        dataTable(1:end,5) = array2table(zeros([height(dataTable) 1]));
    else
        dataTable(:,5) = array2table(append(table2array(dataTable(:,4)), newCompartmentName)); %add COMMON compartment
    end

    dataTable(:,6) = printRxnFormula(updatedModel, table2array(dataTable(:,4)));  %check reaction formula
    dataTable = table2array(dataTable);


% If COMMON mets are provided
else  
    dataTable = OG_mets;  % Add model met names



    dataTable(:,2) = array2table(erase(table2array(dataTable(:,1)), "[e]")); % remove compartment
    dataTable(:,2) = array2table(erase(table2array(dataTable(:,2)), "[c]")); % remove compartment
    dataTable(:,3) = COMMON_mets; %add COMMON compartment
    dataTable(:,4) = Rxns; % get model reaction IDs

    if add_COMMON_EX == 0
        dataTable(1:end,5) = array2table(num2cell(zeros([height(dataTable) 1])));
    else
        dataTable(:,5) = array2table(append(table2array(dataTable(:,4)), newCompartmentName)); %add COMMON compartment
    end

    dataTable(:,6) = printRxnFormula(updatedModel, table2array(dataTable(:,4)));  %check reaction formula
    dataTable = table2array(dataTable);


end
    

%% Change all EX reactions to be met <=> met[i] and add [i] met 
for n = 1:height(dataTable)


    updatedModel = addReaction(updatedModel, convertStringsToChars(dataTable(n,4)),...
        'metaboliteList', {convertStringsToChars(dataTable(n,1)),convertStringsToChars(dataTable(n,3))},...
        'stoichCoeffList', [-1 ; 1], ...
        'lowerBound',-1000, ...
        'upperBound', 1000, ...
        "subSystem", {'Transport'});


disp(n)
        dataTable(n,7) = printRxnFormula(updatedModel, convertStringsToChars(dataTable(n,4)));
        dataTable(n,8) = num2cell(updatedModel.lb(findRxnIDs(updatedModel,convertStringsToChars(dataTable(n,4)))));
        dataTable(n,9) = num2cell(updatedModel.ub(findRxnIDs(updatedModel,convertStringsToChars(dataTable(n,4)))));





%% Execute if addition of met[i] <=> is desired
    if add_COMMON_EX == 1
        updatedModel = addReaction(updatedModel, convertStringsToChars(dataTable(n,5)),...  % Add met[i] <=> []
            'metaboliteList', {convertStringsToChars(dataTable(n,3))},...
            'stoichCoeffList', -1, ...
            'lowerBound', lb(n,1), ...
            'upperBound', ub(n,1), ...
            "subSystem", {'Exchange with the environment'});

        dataTable(n,10) = printRxnFormula(updatedModel, convertStringsToChars(dataTable(n,5)));

        dataTable(n,11) = updatedModel.lb(int64(findRxnIDs(updatedModel,convertStringsToChars(dataTable(n,5)))));
        dataTable(n,12) = updatedModel.ub(int64(findRxnIDs(updatedModel,convertStringsToChars(dataTable(n,5)))));


    end


end
disp(add_COMMON_EX)
dataTable = array2table(dataTable);

%% Rename output table
if add_COMMON_EX == 0
    dataTable.Properties.VariableNames = {'OG_mets' 'No_comp' 'COMMON_mets' 'model.rxns' 'COMMON_rxns' 'OG_rxns_Formulae' 'COMMON_rxns_Formulae' 'COMMON_rxn_lb' 'COMMON_rxn_ub'};
end
if add_COMMON_EX == 1
    dataTable.Properties.VariableNames = {'OG_mets' 'No_comp' 'COMMON_mets' 'model.rxns' 'COMMON_rxns' 'OG_rxns_Formulae' 'COMMON_rxns_Formulae' 'COMMON_rxn_lb' 'COMMON_rxn_ub' 'COMMON_EX_Formulae' 'EX_rxn_lb' 'EX_rxn_ub'};
end





end
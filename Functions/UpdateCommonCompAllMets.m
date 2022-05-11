function [updatedModel, dataTable] = UpdateCommonCompAllMets(model,iCEL_conversionTable, newCompartmentName, add_COMMON_EX, lb, ub, fix_h)
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
 iCEL_conversionTable
 newCompartmentName = '[i]';
 add_COMMON_EX  = 0  ;
 lb = -1000;
 ub = 1000;

 fix_h = 0; 
end


updatedModel = model;






%% Modify all Exchange reactions to be met <=> met[i]

dataTable = iCEL_conversionTable(:, "NEW_EX_names");
for n = 1:height(iCEL_conversionTable)

    updatedModel = addReaction(updatedModel, ...
        char(table2array(iCEL_conversionTable(n,5))) ...
        ,...
        'metaboliteList', {char(table2array(iCEL_conversionTable(n,2))),char(table2array(iCEL_conversionTable(n,3)))},...
        'stoichCoeffList', [ -1 ; 1], ...
        'lowerBound',-1000, ...
        'upperBound', 1000, ...
        "subSystem", {'Transport'});



    dataTable(n,2) = printRxnFormula(updatedModel, table2array(dataTable(n,1)));
    dataTable(n,3) = array2table(updatedModel.lb(findRxnIDs(updatedModel,table2array(dataTable(n,1)))));
    dataTable(n,4) = array2table(updatedModel.ub(findRxnIDs(updatedModel,table2array(dataTable(n,1)))));







    %% Add met[i] <=> [] Exchange reactions

    if add_COMMON_EX == 1

        updatedModel = addReaction(updatedModel, ...
            char(table2array(iCEL_conversionTable(n,"COMMON_EX"))),...
            'metaboliteList', {char(table2array(iCEL_conversionTable(n,"COMMON_mets")))},...
            'stoichCoeffList', -1, ...
            'lowerBound', lb(n,1), ...
            'upperBound', ub(n,1), ...
            "subSystem", {'Exchange with the environment'});

        dataTable(n,5) = printRxnFormula(updatedModel, table2array(iCEL_conversionTable(n,"COMMON_EX")));

        dataTable(n,6) = array2table(updatedModel.lb(findRxnIDs(updatedModel,table2array(iCEL_conversionTable(n,6)))));
        dataTable(n,7) = array2table(updatedModel.ub(findRxnIDs(updatedModel,table2array(iCEL_conversionTable(n,6)))));


    end
end

if add_COMMON_EX == 0
    dataTable.Properties.VariableNames = {'NEW_EX_name' 'NEW_EX_rxnFormula' 'NEW_EX_lb' 'NEW_EX_ub'};
end
if add_COMMON_EX == 1
    dataTable.Properties.VariableNames = {'NEW_EX_name' 'NEW_EX_rxnFormula' 'NEW_EX_lb' 'NEW_EX_ub' 'IEX_rxnFormula' 'IEX_lb' 'IEX_ub'};
end





end
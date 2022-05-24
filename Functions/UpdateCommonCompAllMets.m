function [updatedModel, dataTable] = UpdateCommonCompAllMets(model,CT, newCompartmentName, add_COMMON_EX, lb, ub)
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
 CT
 newCompartmentName = '[i]';
 add_COMMON_EX  = 0  ;
 lb = -1000;
 ub = 1000;

end


updatedModel = model;






%% Modify all Exchange reactions to be met <=> met[i]

dataTable = CT(:, "NEW_IEX_names");
for n = 1:height(CT)

    updatedModel = addReaction(updatedModel, ...
        char(table2array(CT(n,5))) ...
        ,...
        'metaboliteList', {char(table2array(CT(n,2))),char(table2array(CT(n,3)))},...
        'stoichCoeffList', [-1 ; 1], ...
        'lowerBound',-1000, ...
        'upperBound', 1000, ...
        "subSystem", {'Transport'});



    dataTable(n,2) = printRxnFormula(updatedModel, table2array(dataTable(n,1)));
    dataTable(n,3) = array2table(updatedModel.lb(int64(findRxnIDs(updatedModel,table2array(dataTable(n,1))))));
    dataTable(n,4) = array2table(updatedModel.ub(int64(findRxnIDs(updatedModel,table2array(dataTable(n,1))))));







    %% Add met[i] <=> [] Exchange reactions

    if add_COMMON_EX == 1

        updatedModel = addReaction(updatedModel, ...
            char(table2array(CT(n,"NEW_EX_names"))),...
            'metaboliteList', {char(table2array(CT(n,"COMMON_mets")))},...
            'stoichCoeffList', -1, ...
            'lowerBound', lb(n,1), ...
            'upperBound', ub(n,1), ...
            "subSystem", {'Exchange with the environment'});

        dataTable(n,5) = printRxnFormula(updatedModel, table2array(CT(n,"NEW_EX_names")));

        dataTable(n,6) = array2table(updatedModel.lb(findRxnIDs(updatedModel,table2array(CT(n,6)))));
        dataTable(n,7) = array2table(updatedModel.ub(findRxnIDs(updatedModel,table2array(CT(n,6)))));


    end
end

if add_COMMON_EX == 0
    dataTable.Properties.VariableNames = {'NEW_IEX_name' 'NEW_IEX_rxnFormula' 'NEW_IEX_lb' 'NEW_IEX_ub'};
end
if add_COMMON_EX == 1
    dataTable.Properties.VariableNames = {'NEW_IEX_name' 'NEW_IEX_rxnFormula' 'NEW_IEX_lb' 'NEW_IEX_ub' 'EX_rxnFormula' 'EX_lb' 'EX_ub'};
end





end
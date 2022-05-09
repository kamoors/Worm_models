function [Merged_Model, dataTables] = Merge(iCEL_SBML_model,gapseq_model, matching_EX_sheet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MERGE.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USAGE:
%   This functions will take as input the iCEL1314 model and a gapseq model
%   and merge them by creating an external compartment where common EXCHANGE
%   metabolites can be secreted in and taken up from.
%
% [Merged_Model, dataTables] = Merge(iCEL_SBML_model,gapseq_model, matching_EX_sheet)

%% INPUTS:
%    iCEL_SBML_model:     iCEL1314 model from .XML file (Yilmaz et al. 2020)
%    gapseq_model:        gapseq model (Zimmerman et al. 2020)
%    matching_EX_sheet:   STRING containing the name of the .XLSX file 
%                         where the matching information is stored 
%                         (from RStudio script)     

%% OUTPUTS:
%    Merged_Model:     COBRA model structure with gapseq reactions added
%                      to iCEL1314 model
%    dataTables:       Informative data relating to the updated
%                      compartmentents and added reactions /
%                      metabolites. This information is generated mainly
%                      from the activities involving the COMMON
%                      compartment
%
%% AUTHOR:
% Karlis Moors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create ''test'' versions of models
iCEL_test = iCEL_SBML_model;

gapseq_test = gapseq_model;

% Load information table obtained through matching iCEL1314 and gapseq exchange metabolites (in RStudio)
Matching_EX_table = readtable(matching_EX_sheet);

% Find metabolites where the charges match (and
% don't match) between the two models 

ind = table2array(Matching_EX_table(:, 6)) == table2array(Matching_EX_table(:, "iCEL_mets_charge"));

charge_match = Matching_EX_table(ind,:);
charge_Nomatch = Matching_EX_table(ind == 0,:);

charge_match.Properties.VariableNames = {'gapseq_rxns' 'gapseq_rxnNames' 'gapseq_mets' 'gapseq_metFormula' 'gapseq_metNames' 'gapseq_metCharges' 'gapseq_metKEGGID' 'gapseq_metBiGGID' 'combined_match' 'matching_KEGG' 'iCEL_mets_formula' 'iCEL_mets_KEGG' 'iCEL_mets_name' 'iCEL_mets_charge' 'iCEL_EX_rxn' 'iCEL_EX_rxn_names'};
charge_Nomatch.Properties.VariableNames = {'gapseq_rxns' 'gapseq_rxnNames' 'gapseq_mets' 'gapseq_metFormula' 'gapseq_metNames' 'gapseq_metCharges' 'gapseq_metKEGGID' 'gapseq_metBiGGID' 'combined_match' 'matching_KEGG' 'iCEL_mets_formula' 'iCEL_mets_KEGG' 'iCEL_mets_name' 'iCEL_mets_charge' 'iCEL_EX_rxn' 'iCEL_EX_rxn_names'};


% Change charges to charges_noMatch mets
% Create all the variables and the conversion table to account for all necessary entries to deal with metabolites whose charges do no match between iCEL1314 and gapseq models
new_rxn_IDs = array2table(append(table2array(charge_Nomatch(:,"iCEL_EX_rxn")), "_proton"));
new_mets_IDs = array2table(erase(table2array(charge_Nomatch(:,"iCEL_mets_KEGG")), ""));
new_mets_IDs = array2table(erase(table2array(new_mets_IDs), "[e]"));
new_mets_IDs = array2table(erase(table2array(new_mets_IDs), "[c]"));
new_mets_IDs = array2table(append(table2array(new_mets_IDs), "_pr[i]"));
new_COMMON_mets = array2table(erase(table2array(new_mets_IDs), "_pr"));
new_COMMON_rxn_IDs = array2table(append(table2array(charge_Nomatch(:,"iCEL_EX_rxn")), "[i]"));

% Create Coversion Table (to be used in next steps) 
conversion_table = charge_Nomatch(:,"iCEL_mets_KEGG");         % iCEL metabolite IDs (based on matching KEGG identifier)
conversion_table(:,2) = new_mets_IDs;                          % Intermediate metabolite IDs
conversion_table(:,3) = charge_Nomatch(:,"iCEL_mets_name");    % Intermediate metabolite names (same as original)
conversion_table(:,4) = new_rxn_IDs;                           % Rxn IDs for metabolite[e]_pr[i] + h[i] <=> metabolite[e][i]
conversion_table(:,5) =  charge_Nomatch(:,"iCEL_mets_charge"); % Charge in iCEL
conversion_table(:,6) = charge_Nomatch(:,"iCEL_mets_formula"); % Formula in iCEL
conversion_table(:,7) = charge_Nomatch(:,6);                   % Charge in gapseq
conversion_table(:,8) = charge_Nomatch(:,4);                   % Formula in gapseq
conversion_table(:,9) = new_COMMON_mets;                       % metabolite[e]_pr[i] IDs
conversion_table(:,10) = new_COMMON_rxn_IDs;                   % Rxn IDs for metabolite[e][i]  <=> [ ] 
conversion_table(:,11) = charge_Nomatch(:,"iCEL_EX_rxn");      % Original EX rxn IDs
conversion_table(:,12) = charge_Nomatch(:,1);                  % gapseq charge_Nomatch IDs  

conversion_table.Properties.VariableNames = {'iCEL_mets_KEGG' 'new_met_IDs' 'new_met_names' 'new_rxn_IDs' 'iCEL_mets_charge' 'iCEL_mets_formula' 'myb11_metCharges' 'myb11_metFormula' 'new_COMMON_mets' 'new_COMMON_rxn_IDs' 'OG_EX_rxn_IDs' 'myb11_rxns'};

conversion_table = table2array(conversion_table);





% Find the original bounds for the Exchange reactions
IDs_charge_match = findRxnIDs(iCEL_test, table2array(charge_match(:,"iCEL_EX_rxn")));


OG_iCEL_EX = iCEL_test.lb(IDs_charge_match);
OG_iCEL_EX(:,2) = iCEL_test.ub(IDs_charge_match);



OG_gapseq_EX = gapseq_test.lb(findRxnIDs(gapseq_test, table2array(charge_match(:,"gapseq_rxns"))));
OG_gapseq_EX(:,2) = gapseq_test.ub(findRxnIDs(gapseq_test, table2array(charge_match(:,"gapseq_rxns"))));




IDs_charge_noMatch = findRxnIDs(iCEL_test, conversion_table(:,11));


noMatch_iCEL_EX = iCEL_test.lb(IDs_charge_noMatch);
noMatch_iCEL_EX(:,2) = iCEL_test.ub(IDs_charge_noMatch);






% Add all reactions that have matching charges
% Overall aim is to change the original Exchange reaction and add two reactions to the iCEL1314 model
% 1) metabolite[e] <=> metabolite[e]_pr[i]    (Rxn_name EX00000), where EX00000 is the original Exchange reaction name with formula metabolite[e] <=> [ ]
% 2)  metabolite[e]_pr[i] + h[i] <=> metabolite[e][i]           OR          metabolite[e]_pr[i] <=> h[i] + metabolite[e][i]    (Rxn_name EX00000_proton)
% 3)  metabolite[e][i]  <=> [ ]  (Rxn_name EX00000[i])


lb = OG_iCEL_EX(:,1);
ub = OG_iCEL_EX(:,2);

disp("Worm model updating...")



[updatedModel_iCEL, dataTable_iCEL] = UpdateCommonComp(iCEL_test,  charge_match(:, 'iCEL_mets_KEGG'),charge_match(:, 'iCEL_EX_rxn'), "[i]",1,lb, ub);


disp("Worm model updated. dataTable_iCEL generated...")

%  metabolite[e] <=> metabolite[e]_pr[i] + h[i] <=> metabolite[e][i]
lb = noMatch_iCEL_EX(:,1);
ub = noMatch_iCEL_EX(:,2);

disp("metabolite[e] <=> metabolite[e]_pr[i] Protonation reactions being addedto iCEL...")
for i = 2:height(conversion_table)

    
   % add all metabolite[e]_pr[i] mets to updatedModel
updatedModel_iCEL = addMetabolite(updatedModel_iCEL,convertStringsToChars(conversion_table(i,2)), convertStringsToChars(conversion_table(i,3)),convertStringsToChars(conversion_table(i,6)));

   % add all metabolite[e] <=> metabolite[e]_pr[i] reactions to updatedModel
updatedModel_iCEL = addReaction(updatedModel_iCEL, convertStringsToChars(conversion_table(i,11)),...  % Change all EX reactions to be met <=> [] and add [i] met 
            'metaboliteList', {convertStringsToChars(conversion_table(i,1)),convertStringsToChars(conversion_table(i,2))},...
                 'stoichCoeffList', [ -1 ; 1], 'lowerBound', lb(i), 'upperBound', ub(i));

% Check that addition succesful



end

disp("Check Protonation reactions (1:4:end)...")
printRxnFormula(updatedModel_iCEL, {conversion_table(1,11)});


disp("metabolite[e]_pr[i] + h[i] <=> metabolite[e][i] Protonation reactions being addedto iCEL...")
% Add a (de-)protonation reaction to the charge_noMatch metabolites for the iCEL 1314 model
% 
% Add all metabolite[e]_pr[i] + h[i] <=> metabolite[e][i] reactions to updatedModel 
% Add all metabolite[e][i]  <=> [ ]  reactions to updatedModel 

updatedModel_iCEL = UpdateChargeDiff(updatedModel_iCEL, conversion_table,2);


lb_gapseq = OG_gapseq_EX(:,1);
ub_gapseq = OG_gapseq_EX(:,2);

[updatedModel_gapseq, dataTable_gapseq_match] = UpdateCommonComp(gapseq_test, charge_match(:, 'gapseq_mets'), charge_match(:, 'gapseq_rxns'), "[i]", 0, lb_gapseq, ub_gapseq, dataTable_iCEL(:,3));



IDs_charge_noMatch_gapseq = findRxnIDs(gapseq_test, conversion_table(:,12));


noMatch_gapseq_EX = gapseq_test.lb(IDs_charge_noMatch_gapseq);
noMatch_gapseq_EX(:,2) = gapseq_test.ub(IDs_charge_noMatch_gapseq);

[updatedModel_gapseq, dataTable_gapseq_nomatch] = UpdateCommonComp(updatedModel_gapseq, charge_Nomatch(:, 3),  charge_Nomatch(:, 1),"[i]", 0, noMatch_gapseq_EX(:,1), noMatch_gapseq_EX(:,2), array2table(conversion_table(:,9)));


Merged_Model = MergeModelsOnEX(updatedModel_iCEL, updatedModel_gapseq);
Merged_Model.comps = {Merged_Model.comps{1:3} "i" "e0" "c0" "p0"}';
Merged_Model.compNames = {Merged_Model.compNames{1:3} "interaction" "bacterial extracellular space" "bacterial cytosol" "bacterial periplasmic space"}';

Merged_Model.description = char("Merged model iCEL1314 and gapseq");
Merged_Model.modelName = char(strcat(iCEL_SBML_model.modelName,"_",gapseq_model.modelID));
Merged_Model.modelID = char(strcat(iCEL_SBML_model.modelName,"_",gapseq_model.modelID));

dataTables = struct();

dataTables.dataTable_iCEL = dataTable_iCEL;
dataTables.dataTable_gapseq_match = dataTable_gapseq_match;
dataTables.dataTable_gapseq_nomatch = dataTable_gapseq_nomatch;





end
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


%% Find metabolites where the charges match (and don't match) between the two models 



ind = table2array(Matching_EX_table(:, 6)) == table2array(Matching_EX_table(:, "iCEL_mets_charge"));

charge_match = Matching_EX_table(ind,:);
charge_Nomatch = Matching_EX_table(ind == 0,:);


% Add column to charge_match and charge_Nomatch to reflect new EX names

to_append = table2array(charge_match(:,"iCEL_mets_KEGG"));
to_append1 = table2array(charge_match(:,3));

charge_match(:,17) = cellstr(append("IEX_", to_append));
charge_match(:,18) = cellstr(append("EX_", to_append));

charge_match(:,19) = cellstr(append("IEX_", to_append1));
charge_match(:,20) = cellstr(append("EX_", to_append1));

to_append = table2array(charge_Nomatch(:,"iCEL_mets_KEGG"));
to_append1 = table2array(charge_Nomatch(:,3));

charge_Nomatch(:,17) = cellstr(append("IEX_", to_append));
charge_Nomatch(:,18) = cellstr(append("EX_", to_append));

charge_Nomatch(:,19) = cellstr(append("IEX_", to_append1));
charge_Nomatch(:,20) = cellstr(append("EX_", to_append1));


table_names = {'gapseq_rxns' 'gapseq_rxnNames' 'gapseq_mets' 'gapseq_metFormula' 'gapseq_metNames' ...
    'gapseq_metCharges' 'gapseq_metKEGGID' 'gapseq_metBiGGID' 'combined_match' 'matching_KEGG' 'iCEL_mets_formula' ...
    'iCEL_mets_KEGG' 'iCEL_mets_name' 'iCEL_mets_charge' 'iCEL_EX_rxn' 'iCEL_EX_rxn_names' 'NEW_IEX_names' 'NEW_EX_names' ...
    'gapseq_NEW_IEX_names' 'gapseq_NEW_EX_names'};
charge_match.Properties.VariableNames = table_names;
charge_Nomatch.Properties.VariableNames = table_names;


%% Rename Exchange reactions of the iCEL1314 model and prepare further steps
iCEL_test = iCEL_SBML_model;
% Find Exchange reactions
iCEL_Exchange = iCEL_test.rxns(findExcRxns(iCEL_test, 0));
EX_ind = findRxnIDs(iCEL_test, iCEL_Exchange);





% Exclude common EX that have differing charge;
charge_Nomatch_ind = findRxnIDs(iCEL_test, table2array(charge_Nomatch(:,"iCEL_EX_rxn")));
c = setdiff(EX_ind, charge_Nomatch_ind);

% Determine original constraints 
iCEL_lb = iCEL_test.lb(EX_ind);
iCEL_ub = iCEL_test.ub(EX_ind);

% Find associated metabolites
[iCEL_EX_mets,~] = findMetsFromRxns(iCEL_test, iCEL_Exchange);
iCEL_EX_mets = iCEL_EX_mets(:,1);

% Determine metabolite without compartment 
erased_comp = erase(string(iCEL_EX_mets), "[e]");
erased_comp = erase(string(erased_comp), "[c]");
erased_comp = erase(string(erased_comp), "[m]");


% Create conversion table for later use
iCEL_ct = array2table(iCEL_Exchange);
iCEL_ct(:,2) = iCEL_EX_mets;
iCEL_ct(:,3) = cellstr(append(erased_comp,"[i]"));
iCEL_ct(:,4) = cellstr(erased_comp);
iCEL_ct(:,5) = cellstr(append("IEX_", string(iCEL_EX_mets)));
iCEL_ct(:,6) = cellstr(append("EX_", string(iCEL_EX_mets)));

iCEL_ct.Properties.VariableNames = {'OLD_EX_names' 'OG_mets' 'COMMON_mets' 'mets_noComp' 'NEW_IEX_names' 'NEW_EX_names'};




% Rename all Exchange reactions to the form "IEX_metabolite"
iCEL_test.rxns(EX_ind) = table2array(iCEL_ct(:,5));
iCEL_Exchange = iCEL_test.rxns(findExcRxns(iCEL_test, 0));

printRxnFormula(iCEL_test, iCEL_Exchange)


%% Update Common Compartment with iCEL_reactions

[updatedModel_iCEL, dataTable_iCEL_all] = UpdateCommonCompAllMets(iCEL_test,iCEL_ct, '[i]', 1, iCEL_lb, iCEL_ub);

iCEL_Exchange_updated = updatedModel_iCEL.rxns(findExcRxns(updatedModel_iCEL, 0, 0));



[charge_match_mets,~] = findMetsFromRxns(updatedModel_iCEL, table2array(charge_match(:,"NEW_EX_names")));
charge_match_mets = array2table(string(charge_match_mets(:,1)));



[charge_Nomatch_mets,~] = findMetsFromRxns(updatedModel_iCEL, table2array(charge_Nomatch(:,"NEW_EX_names")));
charge_Nomatch_mets = array2table(string(charge_Nomatch_mets(:,1)));
%% Rename Exchange reactions of the gapseq model and prepare further steps
gapseq_test = gapseq_model;
% Find Exchange reactions
gapseq_Exchange = gapseq_test.rxns(findExcRxns(gapseq_test, 0));
EX_ind = findRxnIDs(gapseq_test, gapseq_Exchange);





% Exclude common EX that have differing charge;
charge_Nomatch_ind = findRxnIDs(gapseq_test, table2array(charge_Nomatch(:,"iCEL_EX_rxn")));
c = setdiff(EX_ind, charge_Nomatch_ind);

% Determine original constraints 
gapseq_lb = gapseq_test.lb(EX_ind);
gapseq_ub = gapseq_test.ub(EX_ind);

% Find associated metabolites
[gapseq_EX_mets,~] = findMetsFromRxns(gapseq_test, gapseq_Exchange);
gapseq_EX_mets = gapseq_EX_mets(:,1);

% Determine metabolite without compartment 
erased_comp = erase(string(gapseq_EX_mets), "[e0]");
erased_comp = erase(string(erased_comp), "[c0]");
erased_comp = erase(string(erased_comp), "[p0]");


% Create conversion table for later use
gapseq_ct_all = array2table(gapseq_Exchange);
gapseq_ct_all(:,2) = gapseq_EX_mets;
gapseq_ct_all(:,3) = cellstr(append(erased_comp,"[i]"));
gapseq_ct_all(:,4) = cellstr(erased_comp);
gapseq_ct_all(:,5) = cellstr(append("IEX_", string(gapseq_EX_mets)));
gapseq_ct_all(:,6) = cellstr(append("EX_", string(gapseq_EX_mets)));
gapseq_ct_all(:,7) = array2table(findRxnIDs(gapseq_test, gapseq_Exchange));

gapseq_ct_all.Properties.VariableNames = {'OLD_EX_names' 'OG_mets' 'COMMON_mets' 'mets_noComp' 'NEW_IEX_names' 'NEW_EX_names' 'RxnIDs'};


% Determine metabolite without compartment 
erased_comp = erase(string(table2array(charge_match(:, "gapseq_mets"))), "[e0]");
erased_comp = erase(string(erased_comp), "[c0]");
erased_comp = erase(string(erased_comp), "[p0]");


% Create conversion table with charge_match metabolites / reactions
gapseq_ct_CM = charge_match(:, "gapseq_rxns");
gapseq_ct_CM(:,2) = charge_match(:, "gapseq_mets");
gapseq_ct_CM(:,3) = cellstr(append(erased_comp,"[i]"));
gapseq_ct_CM(:,4) = cellstr(erased_comp);
gapseq_ct_CM(:,5) = cellstr(append("IEX_", table2array(charge_match(:,3))));
gapseq_ct_CM(:,6) = cellstr(append("EX_", table2array(charge_match(:,3))));
gapseq_ct_CM(:,7) = array2table(findRxnIDs(gapseq_test, table2array(gapseq_ct_CM(:,1))));


gapseq_ct_CM.Properties.VariableNames = {'OLD_EX_names' 'OG_mets' 'COMMON_mets' 'mets_noComp' 'NEW_IEX_names' 'NEW_EX_names' 'RxnIDs'};


% Determine metabolite without compartment 
erased_comp = erase(string(table2array(charge_Nomatch(:, "gapseq_mets"))), "[e0]");
erased_comp = erase(string(erased_comp), "[c0]");
erased_comp = erase(string(erased_comp), "[p0]");


% Create conversion table with charge_match metabolites / reactions
gapseq_ct_CNM = charge_Nomatch(:, "gapseq_rxns");
gapseq_ct_CNM(:,2) = charge_Nomatch(:, "gapseq_mets");
gapseq_ct_CNM(:,3) = cellstr(append(erased_comp,"[i]"));
gapseq_ct_CNM(:,4) = cellstr(erased_comp);
gapseq_ct_CNM(:,5) = cellstr(append("IEX_", table2array(charge_Nomatch(:,3))));
gapseq_ct_CNM(:,6) = cellstr(append("EX_", table2array(charge_Nomatch(:,3))));
gapseq_ct_CNM(:,7) = array2table(findRxnIDs(gapseq_test, table2array(gapseq_ct_CNM(:,1))));


gapseq_ct_CNM.Properties.VariableNames = {'OLD_EX_names' 'OG_mets' 'COMMON_mets' 'mets_noComp' 'NEW_IEX_names' 'NEW_EX_names' 'RxnIDs'};


% Determine gapseq metabolites that are not in charge_match and
% charge_Nomatch lists
[diff, IA] = setdiff(gapseq_ct_all(:,7), gapseq_ct_CM(:,7));

ct_all = gapseq_ct_all(IA,:)

[diff, IB] = setdiff(gapseq_ct_all(IA,7), gapseq_ct_CNM(:,7));

ct_all = ct_all(IB,:)

% Rename all Exchange reactions to the form "IEX_metabolite"
gapseq_test.rxns(EX_ind) = table2array(gapseq_ct_all(:,5));
gapseq_Exchange = gapseq_test.rxns(findExcRxns(gapseq_test, 0));

printRxnFormula(gapseq_test, gapseq_Exchange)


%% Update Common Compartment with iCEL_reactions that are not matched

[updatedModel_gapseq, dataTable_gapseq_all] = UpdateCommonCompAllMets(gapseq_test, ct_all, '[i]', 1, gapseq_lb, gapseq_ub);

gapseq_Exchange_updated = updatedModel_gapseq.rxns(findExcRxns(updatedModel_gapseq, 0, 0)); % = 0x1 cell array
% All reactions were changed to met[e0] <-> []  to met[e0] <-> met[i]
% So, there should be no EX reactions in this model

printRxnFormula(updatedModel_gapseq, gapseq_Exchange)


%% Create a conversion data table to add protonation reactions

% Change charges to charges_noMatch mets
% Create all the variables and the conversion table to account for all necessary entries to deal with metabolites whose charges do no match between iCEL1314 and gapseq models
new_rxn_IDs = array2table(append("Proton_", table2array(charge_Nomatch(:,"NEW_EX_names")))); % Create new protonation reaction names
new_mets_IDs = array2table(erase(table2array(charge_Nomatch(:,"iCEL_mets_KEGG")), ""));      % Create new metabolite ids - met_pr[i]
new_mets_IDs = array2table(erase(table2array(new_mets_IDs), "[e]"));
new_mets_IDs = array2table(erase(table2array(new_mets_IDs), "[c]"));
new_mets_IDs = array2table(append(table2array(new_mets_IDs), "_pr[i]"));
new_COMMON_mets = array2table(erase(table2array(new_mets_IDs), "_pr"));
new_COMMON_rxn_IDs = charge_Nomatch(:,"NEW_EX_names");

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
conversion_table(:,11) = charge_Nomatch(:,'NEW_IEX_names');    % Original EX rxn IDs
conversion_table(:,12) = charge_Nomatch(:,'gapseq_NEW_IEX_names');                  % gapseq charge_Nomatch IDs

conversion_table.Properties.VariableNames = {'iCEL_mets_KEGG' 'new_met_IDs' 'new_met_names' 'new_rxn_IDs' 'iCEL_mets_charge' 'iCEL_mets_formula' 'myb11_metCharges' 'myb11_metFormula' 'new_COMMON_mets' 'new_COMMON_rxn_IDs' 'OG_EX_rxn_IDs' 'myb11_rxns'};

conversion_table = table2array(conversion_table);




%% Update Common Compartment with mets that differ in charge

% Find the original bounds for the Exchange reactions
IDs_charge_match = findRxnIDs(iCEL_SBML_model, table2array(charge_match(:,"iCEL_EX_rxn")));


OG_iCEL_EX = iCEL_SBML_model.lb(IDs_charge_match);
OG_iCEL_EX(:,2) = iCEL_SBML_model.ub(IDs_charge_match);


% Find orgiginal bounds for the gapseq model
OG_gapseq_EX = gapseq_model.lb(findRxnIDs(gapseq_model, table2array(charge_match(:,"gapseq_rxns"))));
OG_gapseq_EX(:,2) = gapseq_model.ub(findRxnIDs(gapseq_model, table2array(charge_match(:,"gapseq_rxns"))));



% Fnd original bounds for EX of Nomatch reactions / metabolites
IDs_charge_noMatch = findRxnIDs(iCEL_test, conversion_table(:,11));


noMatch_iCEL_EX = iCEL_SBML_model.lb(int64(IDs_charge_noMatch));
noMatch_iCEL_EX(:,2) = iCEL_SBML_model.ub(IDs_charge_noMatch);





%% Add all reactions that have matching charges

% Overall aim is to change the original Exchange reaction and add two reactions to the iCEL1314 model
% 1) metabolite[e] <=> metabolite[e]_pr[i]    (Rxn_name EX00000), where EX00000 is the original Exchange reaction name with formula metabolite[e] <=> [ ]
% 2)  metabolite[e]_pr[i] + h[i] <=> metabolite[e][i]           OR          metabolite[e]_pr[i] <=> h[i] + metabolite[e][i]    (Rxn_name EX00000_proton)
% 3)  metabolite[e][i]  <=> [ ]  (Rxn_name EX00000[i])


lb = OG_iCEL_EX(:,1);
ub = OG_iCEL_EX(:,2);

disp("Worm model updating...")


% Update the
%[updatedModel_iCEL, dataTable_iCEL] = UpdateCommonComp(updatedModel_iCEL,  charge_match(:, 'iCEL_mets_KEGG'),charge_match(:, 'NEW_IEX_names'), "[i]",1,lb, ub);


disp("Worm model updated. dataTable_iCEL generated...")

%  metabolite[e] <=> metabolite[e]_pr[i] + h[i] <=> metabolite[e][i]
lb = noMatch_iCEL_EX(:,1);
ub = noMatch_iCEL_EX(:,2);

disp("metabolite[e] <=> metabolite[e]_pr[i] Protonation reactions being addedto iCEL...")
for i = 2:height(conversion_table)


    % add all metabolite[e]_pr[i] mets to updatedModel
    updatedModel_iCEL = addMetabolite(updatedModel_iCEL,convertStringsToChars(conversion_table(i,2)), convertStringsToChars(conversion_table(i,3)),convertStringsToChars(conversion_table(i,6)));

    % add all metabolite[e] <=> metabolite[e]_pr[i] reactions to updatedModel
    updatedModel_iCEL = addReaction(updatedModel_iCEL, ...
        convertStringsToChars(conversion_table(i,11)),...  % Change all EX reactions to be met <=> [] and add [i] met
        'metaboliteList', {convertStringsToChars(conversion_table(i,1)),convertStringsToChars(conversion_table(i,2))},...
        'stoichCoeffList', [ -1 ; 1], ...
        'lowerBound', -1000, ...
        'upperBound', 1000);

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

%% Update gapseq model with reactions that are matched

lb_gapseq = OG_gapseq_EX(:,1);
ub_gapseq = OG_gapseq_EX(:,2);

[updatedModel_gapseq, dataTable_gapseq_match] = UpdateCommonComp(updatedModel_gapseq, charge_match(:, 'gapseq_mets'), charge_match(:, 'gapseq_NEW_IEX_names'), "[i]", 0, lb_gapseq, ub_gapseq, charge_match_mets);



IDs_charge_noMatch_gapseq = findRxnIDs(gapseq_model, table2array(charge_Nomatch(:,1)));


lb_noMatch_gapseq = gapseq_model.lb(IDs_charge_noMatch_gapseq);
ub_noMatch_gapseq = gapseq_model.ub(IDs_charge_noMatch_gapseq);

[updatedModel_gapseq, dataTable_gapseq_nomatch] = UpdateCommonComp(updatedModel_gapseq, charge_Nomatch(:, 3),  charge_Nomatch(:, "gapseq_NEW_IEX_names"),"[i]", 0, lb_noMatch_gapseq, ub_noMatch_gapseq, array2table(conversion_table(:,9)));


Merged_Model = MergeModelsOnEX(updatedModel_iCEL, updatedModel_gapseq);
Merged_Model.comps = {Merged_Model.comps{1:3} "i" "e0" "c0" "p0"}';
Merged_Model.compNames = {Merged_Model.compNames{1:3} "interaction" "bacterial extracellular space" "bacterial cytosol" "bacterial periplasmic space"}';

Merged_Model.description = char("Merged model iCEL1314 and gapseq");
Merged_Model.modelName = char(strcat(iCEL_SBML_model.modelName,"_",gapseq_model.modelID));
Merged_Model.modelID = char(strcat(iCEL_SBML_model.modelName,"_",gapseq_model.modelID));


interaction_lb = gapseq_model.lb(findRxnIDs(gapseq_model, table2array(Matching_EX_table(:,1))));

for i = 1:height(charge_match)

Merged_Model = changeRxnBounds(Merged_Model, table2array(charge_match(i, "NEW_EX_names")), interaction_lb(i), 'l');

end





dataTables = struct();

dataTables.dataTable_iCEL = dataTable_iCEL_all;
dataTables.dataTable_gapseq_match = dataTable_gapseq_match;
dataTables.dataTable_gapseq_nomatch = dataTable_gapseq_nomatch;
dataTables.NEW_EX_match = cat(1, charge_match(:,"NEW_EX_names"), charge_Nomatch(:,"NEW_EX_names"));
dataTables.NEW_EX_gapseq = ct_all(:,"NEW_EX_names");
dataTables.NEW_EX_iCEL = iCEL_ct(:,"NEW_EX_names");


end
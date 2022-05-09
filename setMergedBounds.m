function [MergedModel] = setMergedBounds(merged_model,important_EX_on_gapseq ,important_EX_on_iCEL, eff_on_iCEL, eff_on_gapseq, datatable_iCEL, datatable_gapseq, datatable_gapseq_nomatch)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
dataTable_myb11 = datatable_gapseq;
dataTable = datatable_iCEL;
dataTable_myb11_nomatch = datatable_gapseq_nomatch;
eff_on_myb11 = eff_on_gapseq;
MergedModel_test = merged_model;
important_EX_on_myb11 = important_EX_on_gapseq;



merged_EX = MergedModel_test.rxns(findExcRxns(MergedModel_test, 1, 1));



% Get original lb values for iCEL and constrain the merged model COMMON reactions
% with those numbers


lb_from_iCEL = str2double(table2array(dataTable(:,"COMMON_rxn_lb"))) ~= 0;
lb_values_iCEL = str2double(table2array(dataTable(lb_from_iCEL,"COMMON_rxn_lb")));

common_from_iCEL = dataTable(lb_from_iCEL,"COMMON_rxns");


for i = 1:height(common_from_iCEL)
MergedModel_test = changeRxnBounds(MergedModel_test, table2array(common_from_iCEL(i,"COMMON_rxns")), lb_values_iCEL(i), 'l');
end


% Constrain all EX reaction lb of MergedModel_test to 0

MergedModel_test = changeRxnBounds(MergedModel_test, merged_EX, 0, 'l');



MergedModel_test = changeObjective(MergedModel_test, 'BIO0010', 1);



% Get original lb values for myb11 and constrain the merged model COMMON reactions
% with those numbers
lb_from_myb = str2double(table2array(dataTable_myb11(:,"COMMON_rxn_lb"))) ~= 0;
lb_values = str2double(table2array(dataTable_myb11(lb_from_myb,"COMMON_rxn_lb")));
common_from_myb11 = dataTable(lb_from_myb,"COMMON_rxns");

for i = 1:height(common_from_myb11)
MergedModel_test = changeRxnBounds(MergedModel_test, table2array(common_from_myb11(i,"COMMON_rxns")), lb_values(i), 'l');
end

printRxnFormula(MergedModel_test, cellstr(table2array(common_from_myb11)))
printFluxBounds(MergedModel_test, cellstr(table2array(common_from_myb11)))

% Do the same with reactions that needed protonation reactions

lb_from_myb_nomatch = str2double(table2array(dataTable_myb11_nomatch(:,"COMMON_rxn_lb"))) ~= 0;
lb_values_nomatch = str2double(table2array(dataTable_myb11_nomatch(lb_from_myb_nomatch,"COMMON_rxn_lb")));
common_from_myb11_nomatch = dataTable(lb_from_myb_nomatch,"COMMON_rxns");



for i = 1:height(common_from_myb11_nomatch)
MergedModel_test = changeRxnBounds(MergedModel_test, table2array(common_from_myb11_nomatch(i,"COMMON_rxns")), lb_values_nomatch(i), 'l');
end

printRxnFormula(MergedModel_test, cellstr(table2array(common_from_myb11_nomatch)))
printFluxBounds(MergedModel_test, cellstr(table2array(common_from_myb11_nomatch)))




printRxnFormula(MergedModel_test, cellstr(table2array(common_from_iCEL)))
printFluxBounds(MergedModel_test, cellstr(table2array(common_from_iCEL)))

% Change BAC bounds 
MergedModel_test = changeRxnBounds(MergedModel_test, "EXC0050", -0.1, 'l');

% Change BIOMASS EX bounds 
MergedModel_test = changeRxnBounds(MergedModel_test, "EX_cpd11416_c0", -1000, 'l');
MergedModel_test = changeRxnBounds(MergedModel_test, "EX_cpd11416_c0", 1000, 'u');


% Keep the bounds open for reactions that are important for both organisms
printRxnFormula(MergedModel_test, eff_on_myb11)
printFluxBounds(MergedModel_test, eff_on_myb11)
for i = 1:length(important_EX_on_myb11)
MergedModel_test = changeRxnBounds(MergedModel_test, eff_on_myb11(i), important_EX_on_myb11(i), 'l');

end

for i = 1:length(important_EX_on_iCEL)
MergedModel_test = changeRxnBounds(MergedModel_test, eff_on_iCEL(i), important_EX_on_iCEL(i), 'l');

end

MergedModel = MergedModel_test;




end
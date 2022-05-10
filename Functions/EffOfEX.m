function [important_EX_on_myb11 ,important_EX_on_iCEL, eff_on_iCEL, eff_on_myb11] = EffOfEX(merged_model)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

merged_model_1 = merged_model;

merged_EX = merged_model_1.rxns(findExcRxns(merged_model_1, 1, 1));


%Find original Growth rates for both model parts

merged_model_1 = merged_model;
merged_model_1 = changeObjective(merged_model_1, 'EX_cpd11416_c0', 1);
solutions = optimizeCbModel(merged_model_1);
gapseq_growthWT = solutions.f;


merged_model_1 = merged_model;
merged_model_1 = changeObjective(merged_model_1, 'BIO0010', 1);
solutions = optimizeCbModel(merged_model_1);
iCEL_growthWT = solutions.f;

% Find which EX reactions affect the growth of iCEL when myb11 growth,
% 'EX_cpd11416_c0', constrained at gapseq_growthWT
% 

merged_model_1 = merged_model;

merged_model_1 = changeRxnBounds(merged_model_1, 'EX_cpd11416_c0',gapseq_growthWT,"b");


merged_model_1 = changeObjective(merged_model_1, 'BIO0010',1);


[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(merged_model_1, 'FBA',merged_EX);


eff_on_iCEL = delRxn(grRateKO < iCEL_growthWT | isnan(grRateKO));

eff_grR_on_iCEL = grRateKO(grRateKO < iCEL_growthWT | isnan(grRateKO));

printFluxBounds(merged_model_1, eff_on_iCEL);

important_EX_on_iCEL = merged_model_1.lb(findRxnIDs(merged_model_1,eff_on_iCEL));

% Find which EX reactions affect the growth of myb11 when iCEL growth,
% 'BIO0100', constrained at 0.0692
% 

merged_model_1 = merged_model;

merged_model_1 = changeRxnBounds(merged_model_1, 'BIO0010',iCEL_growthWT,"b");


merged_model_1 = changeObjective(merged_model_1, 'EX_cpd11416_c0',1);


[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(merged_model_1, 'FBA',merged_EX);


eff_on_myb11 = delRxn(grRateKO < gapseq_growthWT | isnan(grRateKO));

eff_grR_on_myb11 = grRateKO(grRateKO < gapseq_growthWT | isnan(grRateKO));

printFluxBounds(merged_model_1, eff_on_myb11);



important_EX_on_myb11  = merged_model_1.lb(findRxnIDs(merged_model_1, eff_on_myb11)); 




end
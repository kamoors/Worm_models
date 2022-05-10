function [output] = FVA_custom(Merged_model, i)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

changeCobraSolver("gurobi")
merged_EX = Merged_model.rxns(findExcRxns(Merged_model, 1, 1));
FVA_minFLUX = table();
FVA_maxFLUX = table();
del_rxn = table();

output = struct();

iCEL_myb11_FVA = Merged_model;

iCEL_myb11_FVA = changeRxnBounds(iCEL_myb11_FVA, 'EX_cpd11416_c0',1.644,"b"); % Constrain the gapseq biomass reaction to optimal value


iCEL_myb11_FVA = changeRxnBounds(iCEL_myb11_FVA, iCEL_myb11_FVA.rxns(i),0,"b"); % Constrain each reaction in the model


[minFlux, maxFlux] = fluxVariability(iCEL_myb11_FVA, 50,'rxnNameList', merged_EX);

FVA_minFLUX = array2table(minFlux);
FVA_maxFLUX = array2table(maxFlux);
del_rxn = iCEL_myb11_FVA.rxns(i);

output.FVA_minFLUX = FVA_minFLUX;
output.FVA_maxFLUX = FVA_maxFLUX;
output.del_rxn = del_rxn;

i

end
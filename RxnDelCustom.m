function [output_Table] = RxnDelCustom(model,obj_list, rxn_list)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
arguments
model
obj_list
rxn_list = model.rxns


end




RXNDEL = table();

output = cell([length(obj_list) 3]);


for i = 1:length(obj_list)
tic

model = changeObjective(model, obj_list(i),1);

[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model, 'FBA',rxn_list);



RXNDEL(:,1) = table(grRatio);
RXNDEL(:,2) = table(grRateKO);
RXNDEL(:,3) = table(hasEffect);
RXNDEL(:,4) = table(delRxn);



RXNDEL.Properties.VariableNames = {'grRatio' 'grRateKO'  'hasEffect' 'delRxn'};

output{i,1} = RXNDEL;
output{i,2} = grRateWT;
output{i,3} = fluxSolution;


toc
end


output_Table = output;

end
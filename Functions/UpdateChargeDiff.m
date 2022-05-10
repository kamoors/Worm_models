function [updatedModel_iCEL] = UpdateChargeDiff(updatedModel_iCEL,conversion_table, starting_row)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

arguments

updatedModel_iCEL
conversion_table
starting_row = 1

end


disp("Starting the update for reactions with differing metabolite charges...")

k = -1;

for i = starting_row:height(conversion_table)

printRxnFormula(updatedModel_iCEL, {conversion_table(i,4)})
disp(strcat("iCEL charge: ",conversion_table(i,5)))
disp(strcat("Gapseq charge: ",conversion_table(i,7)))


syms x


S = solve(conversion_table(i,5) + x == conversion_table(i,7),x);




% Add all metabolite[e]_pr[COMMON] + h[COMMON] <=> metabolite[e][COMMON] reactions to updatedModel 

updatedModel_iCEL = addReaction(updatedModel_iCEL, ...
    convertStringsToChars(conversion_table(i,4)),...  
            'metaboliteList', {convertStringsToChars(conversion_table(i,2)),convertStringsToChars("h[i]"), convertStringsToChars(conversion_table(i,9))},...
            'stoichCoeffList', [ -1 ; int64(k*S); 1], ...
            'lowerBound', -1000, ...
            'upperBound', 1000, ...
            "subSystem", {'Protonation'});


% Add all metabolite[e][COMMON]  <=> [ ]  reactions to updatedModel
updatedModel_iCEL = addReaction(updatedModel_iCEL, ...
    convertStringsToChars(conversion_table(i,10)),...   
            'metaboliteList', {convertStringsToChars(conversion_table(i,9))},...
                 'stoichCoeffList', -1, ...
                 'lowerBound', 0, ...
                 'upperBound', 1000, ...
                 "subSystem", {'Exchange with the environment'});


printRxnFormula(updatedModel_iCEL, {conversion_table(i,4)})

disp("New Stoichiometry:")
disp([ -1 ; int64(k*S); 1])


end


end
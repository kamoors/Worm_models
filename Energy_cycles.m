function [out_table] = Energy_cycles(MergedModel,min_reactions, previous_reactions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

arguments

MergedModel
min_reactions
previous_reactions = table()
end
merged_EX = MergedModel.rxns(findExcRxns(MergedModel, 1, 1));  % Find merged model EX reactions


% Initialize results table
out_table = cell([height(min_reactions) 1]);

for i = 1:height(min_reactions)

    % Constrain all merged_EX reactions to 0
iCEL_myb11_0 = changeRxnBounds(MergedModel, merged_EX,0,"l");

iCEL_myb11_0 = changeObjective(iCEL_myb11_0,'RCC0005',1);

% Change ATPM reaction upper bound to ''inf''
iCEL_myb11_0 = changeRxnBounds(iCEL_myb11_0, 'RCC0005',1000000,"u");

if ~isempty(previous_reactions)
    disp("previous_reactions is not empty")
% Constrain the previously checked reactions to 0
iCEL_myb11_0 = changeRxnBounds(iCEL_myb11_0, table2array(previous_reactions),0,"b");
printFluxBounds(iCEL_myb11_0, table2array(previous_reactions))
end

% Constrain each reaction from the min_reactions list to 0
% Then do RDA and check the next list
iCEL_myb11_0 = changeRxnBounds(iCEL_myb11_0, table2array(min_reactions(i,4)),0,"b");
disp("This reaction has been constrained:")
table2array(min_reactions(i,4))

disp("Reactions and their bounds in the min_reactions list:")
printFluxBounds(iCEL_myb11_0, table2array(min_reactions(:,4))) 




out = RxnDelCustom(iCEL_myb11_0, {'RCC0005'});

out_table{i,1} = out{1,1}; % Reaction Deletion output: grRatio | grRateKO | hasEffect | delRxn 





reactions = out_table{i,1};




% Get the deleted reactions that produce the lowest ATPM flux
ind = table2array(reactions(:, "delRxn")) ~= "RCC0005";
reactions = reactions(ind,:);
min_flux = min(round(table2array(reactions(:,2)), 3));

% Reactions that produce the lowest ATPM flux in one table
min_reactions_2 = reactions(round(table2array(reactions(:,2)),3) == min_flux, :);



% Reactions that produce the lowest ATPM flux in one table
out_table{i,2} = min_flux; % Minimal achieved ATPM flux
out_table{i,3} = table2array(min_reactions(i,4));  % Reaction being deleted
out_table{i,4} = min_reactions_2; % New reaction array that cause min flux ATPM


end
    

end

function [matchedFromName,matchedFromFormula] = MatchEXrxns(gapseq_model, other_model)
%MatchEXrxns --- This function takes a GapSeq model and one other GSMM, and
%attempts to match the Exchange reaction metabolites based on name and
%formula


%INPUT:
% gapseq_model = Gapseq model 
% other_model = other GSMM (e.g. C. elegans model iCEL1314 or Wormjam)



%OUTPUT:

% matchedFromName ===> List of metabolites that have a direct name match
% with the metNames from the other mode. 
%%%%%%WARNING: this function removes the compartment suffix from the gapseq
%model metNames!


%matchedFromName columns - [gapseq.metNames, 0 or 1 (can ignore),
%gapseq.formula, other.metNames]

%matchedFromFormula ===> 

%matchedFromFormula columns - [gapseq.metFormula, gapseq.metNames,
%other.metNames based on the matches between the formulae]



%%
[gapseq_selExc, gapseq_selUpt] = findExcRxns(gapseq_model);

[other_selExc, other_selUpt] = findExcRxns(other_model);

x = int64( ...
    findMetIDs(other_model, ...
    findMetsFromRxns(other_model,other_model.rxns(other_selExc))));

EX_mets_2 = other_model.metNames(x);
EX_mets_2(:,2) = num2cell(zeros([length(EX_mets_2) 1]));
EX_mets_2(:,3) = other_model.metFormulas(x);



y = int64( ...
    findMetIDs(gapseq_model, ...
    findMetsFromRxns(gapseq_model,gapseq_model.rxns(gapseq_selExc))));


EX_mets_gapseq = gapseq_model.metNames(y);
EX_mets_gapseq(:,2) = num2cell(zeros([201 1]));
EX_mets_gapseq(:,3) = gapseq_model.metFormulas(y);


clear x y
%%
% Fuzzy search in metNames between iCEL and myb11
%Create new array for myb11 EX that does not contain suffix in the name
newStr = erase(EX_mets_gapseq(:,1),"-e0");
Ex_matching = string([201,4]);
for n = 1:length(newStr)
    
   [d,A] =  fzsearch(EX_mets_2(:,1), newStr(n));

   
 % Length of d will be 2 when there is an exact match 


if length(d) == 2
    Ex_matching(n,1) = newStr(n);
    Ex_matching(n,2) = EX_mets_2(d(1,2));
    
    EX_mets_gapseq(n,2) = num2cell(1);
    EX_mets_2(d(1,2),2) = num2cell(1);
    
     EX_mets_gapseq(n,4) = EX_mets_2(d(1,2),1);
    
else
   EX_mets_gapseq(n,2) = num2cell(0); 

end

    

   
   
end

disp('Number of matched EX reactions:')
disp('GapseqModel:')
sum(cell2mat(EX_mets_gapseq(:,2)), 'All')


disp('Other Model:')
sum(cell2mat(EX_mets_2(:,2)), 'All')


x = cell2mat(EX_mets_gapseq(:,2)) == 0;
x1 = cell2mat(EX_mets_gapseq(:,2)) == 1;
matchedFromName = EX_mets_gapseq(x1,:);
narrowGapseq_EX = EX_mets_gapseq(x,:);

y = cell2mat(EX_mets_2(:,2)) == 0;
narrow2_EX = EX_mets_2(y,:);

clear x y d A

%%
newStr = narrowGapseq_EX(:,3);
Ex_matching = cell(123);
for n = 1:length(newStr)
    
   d = contains(narrow2_EX(:,3), newStr(n));

   



if sum(d, 'all') > 0
    Ex_matching{n,1} = string(newStr(n));
    Ex_matching{n,2} = narrowGapseq_EX(n,1);
    Ex_matching{n,3} = narrow2_EX(d,1);
    
    narrowGapseq_EX(n,2) = num2cell(1);
    %narrowiCEL_EX(d(1,2),2) = num2cell(1);

end
   
end


matchedFromFormula = Ex_matching;
end


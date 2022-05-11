function [results] = CheckFutileCycles(merged_model,iCEL_model,gapseq_model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C. elegans biomass reactions
% 'BIO0100'	Body (including storage)
% 'BIO0101'	Body (excluding storage)
% 'BIO0102'	Embryo (biosynthesis in germline)
% 'BIO0103'	Embryo (growth in eggs)
% 'BIO0010'	Any mixture

%% C. elegans energetics
% 'RCC0005'	Non-growth associated maintenance



%% gapseq biomass reaction (sink):
% 'EX_cpd11416_c0'

% gapseq energetics

% "rxn00062" ATP maintenance requirement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results = table();

%% Find EX reactions 
merged_EX = merged_model.rxns(findExcRxns(merged_model, 1, 1));
iCEL_EX = iCEL_model.rxns(findExcRxns(iCEL_model, 1, 1));
gapseq_EX = gapseq_model.rxns(findExcRxns(gapseq_model, 1, 1));
%%%%%%%%%%%%%%%%
%% iCEL1314 model
disp("(1) iCEL1314 original biomass ('BIO0010' Any mixture):")
solutions = optimizeCbModel(iCEL_model);
iCEL_OG_biomass = solutions.f

results(1,1) = array2table(solutions.f);

disp("(2) iCEL1314 biomass; all EX constrained to 0:")
iCEL_model = changeRxnBounds(iCEL_model, iCEL_EX,0,"b");
solutions = optimizeCbModel(iCEL_model);
solutions.f
results(1,2) = array2table(solutions.f);

disp("(3) iCEL1314 NGAM; all EX constrained to 0:")
iCEL_model = changeObjective(iCEL_model, 'RCC0005', 1);
solutions = optimizeCbModel(iCEL_model);
solutions.f
results(1,3) = array2table(solutions.f);

%%%%%%%%%%%%%%
%% gapseq model
disp(strcat("(4) ",gapseq_model.modelID, " original biomass (EX_cpd11416_c0):"))
solutions = optimizeCbModel(gapseq_model);
gapseq_OG_bio = solutions.f
results(1,4) = array2table(solutions.f);


disp(strcat("(5) ",gapseq_model.modelID, " biomass; all EX constrained to 0:"))
gapseq_model = changeRxnBounds(gapseq_model, gapseq_EX,0,"b");
solutions = optimizeCbModel(gapseq_model);
solutions.f
results(1,5) = array2table(solutions.f);

disp(strcat("(6) ",gapseq_model.modelID, " NGAM; all EX constrained to 0:"))
gapseq_model = changeObjective(gapseq_model, 'rxn00062_c0', 1);
gapseq_model = changeRxnBounds(gapseq_model, gapseq_EX,0,"b");
solutions = optimizeCbModel(gapseq_model);
solutions.f
results(1,6) = array2table(solutions.f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merged model - model-part check

%% With iCEL biomass


disp("(7) Merged model original iCEL biomass ('BIO0010' Any mixture):")
merged_model_1 = merged_model;
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,7) = array2table(solutions.f);


disp("(8) Merged model iCEL biomass ('BIO0010' Any mixture); iCEL EX constrained to 0:")
merged_model_1 = changeRxnBounds(merged_model_1, iCEL_EX,0,"b");
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,8) = array2table(solutions.f);

disp(strcat("(9) Merged model iCEL biomass ('BIO0010' Any mixture); ", gapseq_model.modelID, " EX constrained to 0:"))
merged_model_1 = merged_model;
merged_model_1 = changeRxnBounds(merged_model_1, gapseq_EX,0,"b");
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,9) = array2table(solutions.f);

disp("(10) Merged model iCEL biomass ('BIO0010' Any mixture); Merged model EX constrained to 0:")
merged_model_1 = merged_model;
merged_model_1 = changeRxnBounds(merged_model_1, merged_EX,0,"b");
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,10) = array2table(solutions.f);

disp("(11) Merged model iCEL NGAM (RCC0005); Merged model EX constrained to 0:")
merged_model_1 = merged_model;
merged_model_1 = changeObjective(merged_model_1, 'RCC0005', 1);
merged_model_1 = changeRxnBounds(merged_model_1, merged_EX,0,"b");
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,11) = array2table(solutions.f);


%% With gapseq biomass
disp(strcat("(12) Merged model original ", gapseq_model.modelID, " biomass (EX_cpd11416_c0)"))
merged_model_1 = merged_model;
merged_model_1 = changeRxnBounds(merged_model_1, "EXC0050",-0.1,"l");
merged_model_1 = changeRxnBounds(merged_model_1, "EX_cpd11416_c0",1000,"u");
merged_model_1 = changeRxnBounds(merged_model_1, "EX_cpd11416_c0",-1000,"l");
merged_model_1 = changeObjective(merged_model_1, 'EX_cpd11416_c0', 1);

solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,12) = array2table(solutions.f);


disp(strcat("(13) Merged model ", gapseq_model.modelID, " biomass (EX_cpd11416_c0);", gapseq_model.modelID, " EX constrained to 0:"))
merged_model_1 = merged_model;
merged_model_1 = changeObjective(merged_model_1, 'EX_cpd11416_c0', 1);
merged_model_1 = changeRxnBounds(merged_model_1, gapseq_EX,0,"b");
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,13) = array2table(solutions.f);

disp(strcat("(14) Merged model ", gapseq_model.modelID, " biomass (EX_cpd11416_c0); iCEL1314 EX constrained to 0:"))
merged_model_1 = merged_model;
merged_model_1 = changeObjective(merged_model_1, 'EX_cpd11416_c0', 1);
merged_model_1 = changeRxnBounds(merged_model_1, iCEL_EX,0,"b");
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,14) = array2table(solutions.f);

disp(strcat("(15) Merged model ", gapseq_model.modelID, " biomass (EX_cpd11416_c0); merged EX constrained to 0:"))
merged_model_1 = merged_model;
merged_model_1 = changeObjective(merged_model_1, 'EX_cpd11416_c0', 1);
merged_model_1 = changeRxnBounds(merged_model_1, merged_EX,0,"b");
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,15) = array2table(solutions.f);

disp(strcat("(16) Merged model ", gapseq_model.modelID, " NGAM (rxn00062_c0); merged EX constrained to 0:"))
merged_model_1 = merged_model;
merged_model_1 = changeObjective(merged_model_1, 'rxn00062_c0', 1);
merged_model_1 = changeRxnBounds(merged_model_1, merged_EX,0,"b");
solutions = optimizeCbModel(merged_model_1);
solutions.f
results(1,16) = array2table(solutions.f);


results.Properties.RowNames = strcat(gapseq_model.modelID,"");

results.Properties.VariableNames = {'iCEL_biomass' 'iCEL_biomass_EX_0' 'iCEL_NGAM_EX_0' 'Gapseq_biomass' 'Gapseq_biomass_EX_0' 'Gapseq_NGAM_EX_0' 'Merged_iCELbiomass' 'Merged_iCELbiomass_iCEL_EX_0' 'Merged_iCELbiomass_gapseq_EX_0' 'Merged_iCELbiomass_Merged_EX_0' 'Merged_iCELNGAM_Merged_EX_0' 'Merged_gapseqbiomass' 'Merged_gapseqbiomass_gapseq_EX_0' 'Merged_gapseqbiomassiCEL_EX_0' 'Merged_gapseqbiomass_Merged_EX_0' 'Merged_gapseqNGAM_Merged_EX_0'}
end



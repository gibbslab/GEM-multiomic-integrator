close all force; clear variables; clc
%% Startup the COBRA Toolbox
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver ('gurobi', 'all', 1); % For large models

%% Import the proteome and the model
abundance = readtable('sup_material_8.xlsx');
load('Recon3D_301.mat')

%% 1. Get the list of ProteinECNumbers
if isfile('ProteinECNumbers.mat')
    load('ProteinECNumbers.mat', 'ProteinECNumbers');
else
    ProteinECNumbers = getProteinCodes(abundance);
    save('ProteinECNumbers.mat', 'ProteinECNumbers');
end

%% 2. Take rxnECNumbers and format correctly 
Recon3D.rxnECNumbers = regexprep(Recon3D.rxnECNumbers,";$","");
Recon3D.rxnECNumbers = strrep(Recon3D.rxnECNumbers,";"," & ");
Recon3D.rxnECNumbers = strrep(Recon3D.rxnECNumbers,","," & ");
Recon3D.rxnECNumbers = strrep(Recon3D.rxnECNumbers,"or","|");
Recon3D.rxnECNumbers = strrep(Recon3D.rxnECNumbers,"EC:","");
Recon3D.rxnECNumbers = strrep(Recon3D.rxnECNumbers,"TCDB:","");


%% Find out how many reactions have a EC Numer asociated
N_filled = find(~cellfun(@isempty,Recon3D.rxnECNumbers));
percent_of_filled = length(N_filled)/length(Recon3D.rxnECNumbers);

Recon3D = makePrules(Recon3D);

%% Map reactions abundances
FinalAbundances = table();
abundance = addvars(abundance, ProteinECNumbers, 'After', 1);
for k=3:length(abundance.Properties.VariableNames)
    abundanceToMap =  abundance(:, [2 k]);
    abundanceToMap.Properties.VariableNames = {'id' 'value'};
    
    % Remove missing values
    abundanceToMap = rmmissing(abundanceToMap);
    abundanceToMap = abundanceToMap(~cellfun('isempty', abundanceToMap.id), :);
    
    [abundanceRxns parsedPR] = mapAbundanceToReactions(Recon3D, abundanceToMap);
    FinalAbundances = addvars(FinalAbundances, abundanceRxns);
end

expresion = readtable('DESeq2_normalised_counts.xls');
%% 1. Convert Ensembl to Entrez with https://www.ensembl.org/biomart/martview
f = fopen('mart_export.txt');
C = textscan(f, '%s %s', 'HeaderLines', 1);
fclose(f);
EnsemblIDs = regexprep(expresion.(1),"\.[0-9]*","");
EnrezIDs = cell(numel(EnsemblIDs), 1);
for i=1:numel(EnsemblIDs)
    match = strcmp(C{1},EnsemblIDs(i));
    if sum(match) > 0
        EnrezIDs{i} = C{2}{match};
    end
end
expresion = addvars(expresion, EnrezIDs, 'After', 1);

%% Map expresions for each reaction
Recon3D.genes = regexprep(Recon3D.genes,"\.[0-9]*","");

expresion_columns = expresion.Properties.VariableNames;
expressionRxns1 = cell(length(expresion_columns) - 2, 1);
FinalExpresion = table();
for k=3:length(expresion_columns)
    expressionToMap =  expresion(:, [2 k]);
    expressionToMap.Properties.VariableNames = {'gene' 'value'};
    expressionToMap.gene(cellfun('isempty', expressionToMap.gene)) = {' '};
    [expressionRxns parsedGPR] = mapExpressionToReactions(Recon3D, expressionToMap);
    FinalExpresion = addvars(FinalExpresion, expressionRxns);
end
FinalExpresion.Properties.VariableNames = expresion.Properties.VariableNames(3:end);
FinalAbundances.Properties.VariableNames = abundance.Properties.VariableNames(3:end);

%% Dimentionality reduction with PCA
X = [FinalAbundances,FinalExpresion];
Y = table2array(X);
Y(Y == -1) = NaN;
out = normalize(Y(any(~isnan(Y), 2), :));
[coeff,pca_scores,latent,tsquared,explained,mu] = pca(out);
reducedDimension = coeff(:,1:2);

out(isnan(out)) = 0;
reducedFeatureMatrix = (out * reducedDimension) * latent(1:2) / sum(latent(1:2));
red = ones(size(X,1), 1) * -1;
red(any(table2array(X) ~= -1, 2)) = reducedFeatureMatrix;

%% Model reduction
options.solver = 'INIT';
options.weights = red;
options.runtime = 14400;

astrociteModel = createTissueSpecificModel(Recon3D, options);
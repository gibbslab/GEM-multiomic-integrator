close all force; clear variables; clc
%% Startup the COBRA Toolbox
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver('glpk', 'all');

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

%% 2. Extract expresions for each reaction

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
FinalExpresion1 = FinalExpresion;
FinalAbundances1 = FinalAbundances;

% 
% X = [FinalAbundances,FinalExpresion];%X = [FinalAbundances,FinalExpresion(:,ExIndex)];
% Y = table2array(X);
% Y(Y == -1) = NaN;
% out = Y(any(~isnan(Y), 2), :);
% [coeff,pca_scores,latent,tsquared,explained,mu] = pca(out);
% 
% ExIndex = [8 7 17 18 27 28 6 5 15 16 25 26 9 10 19 20 29 30];
% scores = {};
% combinedValues = zeros(numel(Recon3D.rxns), numel(ExIndex));
% for i=1:numel(FinalAbundances.Properties.VariableNames)
%     expr = table2array(FinalExpresion(:,ExIndex(i)));
%     abun = table2array(FinalAbundances(:, i));
%     values = [expr, abun];
%     X = [FinalAbundances,FinalExpresion(:,ExIndex)]
%     [coeff,pca_scores,latent,tsquared,explained,mu] = pca(table2array(X));
%     scores{i} = pca_scores;
% 
% end
% 
% 
% x = [-0.0678   -0.6460    0.5673    0.5062;
%    -0.6785   -0.0200   -0.5440    0.4933;
%     0.0290    0.7553    0.4036    0.5156;
%     0.7309   -0.1085   -0.4684    0.4844];
% [coeff,pca_scores,latent,tsquared,explained,mu] = pca(normalize(x));
% 
% 
% 
% options.solver = 'iMAT';
% options.expressionRxns = combinedValues(1, :);
% options.threshold_lb = -1000;
% options.threshold_ub = 1000;
% 
% astrociteModel = createTissueSpecificModel(Recon3D, options);
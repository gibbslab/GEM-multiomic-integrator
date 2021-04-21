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
FinalAbundances.Properties.VariableNames = abundance.Properties.VariableNames(2:end);

exVal = FinalExpresion{40,[6 7 16 17 26 27]};
abVal = FinalAbundances(44,[8 9 10 11 12 13]);
abVal =  [27.4479 27.3967 27.2953 27.8241 27.5015 27.5533]
scatter(exVal([1 2]), abVal([1 2]))
hold on
scatter(exVal([3 4]), abVal([3 4]))
scatter(exVal([5 6]), abVal([5 6]))

close all force; clear variables; clc
%% Startup the COBRA Toolbox
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver('glpk', 'all');
%% Load the model
model=readCbModel('Astrocyte_Osorio2019.xml');
abundance = readtable('abundance.xlsx');

%% Get the list of ProteinECNumbers
if isfile('ProteinECNumbers.mat')
    load('ProteinECNumbers.mat', 'ProteinECNumbers');
else
    ProteinECNumbers = getProteinCodes(abundance);
    save('ProteinECNumbers.mat', 'ProteinECNumbers');
end


% Go to https://www.uniprot.org/uploadlists/ select Entrez to UniprotKB
% Enable the EC numbers column and download as excel
ECnumbers_osorio=readtable('ECnumbers_Osorio.xlsx');
ECNumbers=table2array(ECnumbers_osorio(:,'ECNumber'));

ECNumbers = regexprep(ECNumbers,';$','');
ECNumbers = strrep(ECNumbers,';',' & ');
ECNumbers = strrep(ECNumbers,',',' & ');
ECNumbers = strrep(ECNumbers,'or','|');
ECNumbers = strrep(ECNumbers,'EC:','');

model = makePrules(model, ECNumbers);

%% Map reactions abundances
FinalAbundances = table();
abundance = addvars(abundance, ProteinECNumbers, 'After', 1);
for k=3:length(abundance.Properties.VariableNames)
    abundanceToMap =  abundance(:, [2 k]);
    abundanceToMap.Properties.VariableNames = {'id' 'value'};
    
    % Remove missing values
    abundanceToMap = rmmissing(abundanceToMap);
    abundanceToMap = abundanceToMap(~cellfun('isempty', abundanceToMap.id), :);
    
    [abundanceRxns parsedPR] = mapAbundanceToReactions(model, abundanceToMap);
    FinalAbundances = addvars(FinalAbundances, abundanceRxns);
end

expresion = readtable('expression.xls');
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

model.genes = regexprep(model.genes,"\.[0-9]*","");

expresion_columns = expresion.Properties.VariableNames;
expressionRxns1 = cell(length(expresion_columns) - 2, 1);
FinalExpresion = table();
for k=3:length(expresion_columns)
    expressionToMap =  expresion(:, [2 k]);
    expressionToMap.Properties.VariableNames = {'gene' 'value'};
    expressionToMap.gene(cellfun('isempty', expressionToMap.gene)) = {' '};
    [expressionRxns parsedGPR] = mapExpressionToReactions(model, expressionToMap);
    FinalExpresion = addvars(FinalExpresion, expressionRxns);
end
FinalExpresion.Properties.VariableNames = expresion.Properties.VariableNames(3:end);
FinalAbundances.Properties.VariableNames = abundance.Properties.VariableNames(3:end);
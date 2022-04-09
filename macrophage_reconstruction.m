close all force; clear variables; clc
%% Startup the COBRA Toolbox
addpath('cobratoolbox')
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver ('gurobi', 'all', 1); % For large models
clc
addpath('functions')

%% Import the proteome and the model
abundance = readtable('data/Proteins-processed.csv');
expresion = readtable('data/expression_data.csv');

%% Get the list of ProteinECNumbers
proteinECNumbers_FileName = 'data/ProteinECNumbers.mat';
if isfile(proteinECNumbers_FileName)
    load(proteinECNumbers_FileName, 'ProteinECNumbers');
else
    ProteinECNumbers = getProteinCodes(abundance);
    save(proteinECNumbers_FileName, 'ProteinECNumbers');
end

%% Map omic data to Recon3D reactions
load('Recon3D_301.mat')
Recon3D = harmonizeReconECnumbers(Recon3D);
Recon3D = makePrules(Recon3D);


% Find out how many reactions have a EC Numer asociated
N_filled = find(~cellfun(@isempty,Recon3D.rxnECNumbers));
percent_of_filled = length(N_filled)/length(Recon3D.rxnECNumbers);


%% Map omic data to Recon3D reactions
disp('Mapping omic data to Recon3D reactions ...')

abundanceTable_FileName = 'out/abundanceTable.mat';
expressionTable_FileName = 'out/expressionTable.mat';
if isfile(expressionTable_FileName) && isfile(abundanceTable_FileName)
    load(abundanceTable_FileName, 'abundanceTable');
    load(expressionTable_FileName, 'expressionTable');
else
    profile on
    abundanceTable = mapAbundance(Recon3D, abundance, ProteinECNumbers);
    expresion = EnsemblToEntrez(expresion, 'data/mart_export.txt');
    expressionTable = mapExpression(Recon3D, expresion);
    profile viewer

    save(abundanceTable_FileName, 'abundanceTable');
    save(expressionTable_FileName, 'expressionTable');
end

%% Transcriptomic and proteomic data integration
disp('Integrating Transcriptomic and Proteomic data ...')
omicIntegratedData = omicIntegrationPCA(abundanceTable, expressionTable);
% omicIntegratedData = log(abs(min(omicIntegratedData)) + omicIntegratedData);

%% Model reconstruction
disp('Creating macrophage Specific Model ...')
options.solver = 'iMAT';
options.expressionRxns = red;
options.threshold_lb = 0.5;
options.threshold_ub = -1;

macrophageModel = createTissueSpecificModel(Recon3D, options);
macrophageModel.description = 'MendozamacrophageModel';
macrophageModel.rxnECNumbers = {macrophageModel.rxnECNumbers{:}}';
save('out/dirtymacrophageModel.mat', 'macrophageModel');

%% Contextualization
[~,idx] = ismember(macrophageModel.rxns, Recon3D.rxns);
omicIntegratedDataMacrophage = omicIntegratedData(idx);
omicIntegratedDataMacrophage(omicIntegratedDataMacrophage == -1) = 0;

model =  exp2flux(macrophageModel, omicIntegratedDataMacrophage);

FBAsolution = optimizeCbModel(model ,'max');
fprintf("The FBA flux after constrain to DMEM %f\n", FBAsolution.f)
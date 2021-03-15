clear variables; clc
%% Startup the COBRA Toolbox
addpath('/home/nico/Documentos/UNAL/Tesis/COBRA/cobratoolbox','-end')
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver('glpk', 'all');

%% Import the proteome and the model
abundance = readtable('sup_material_9.xlsx');
load('Recon3D_301.mat')

%% Get the list of ProteinECNumbers
ProteinECNumbers = mapAbundanceToReactions(Recon3D, abundance);
save('ProteinECNumbers.mat', 'ProteinECNumbers')
% load('ProteinECNumbers.mat')

% How matches in simple cell lists
find(strcmp(Recon3D.rxnECNumbers, '1.4.3.6'));

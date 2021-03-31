close all force; clear variables; clc
%% Startup the COBRA Toolbox
addpath('/home/nico/Documentos/UNAL/Tesis/COBRA/cobratoolbox','-end')
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver('glpk', 'all');

%% Import the proteome and the model
abundance = readtable('sup_material_9.xlsx');
load('Recon3D_301.mat')

%% Find out how many reactions have a EC Numer asociated
N_filled = find(~cellfun(@isempty,Recon3D.rxnECNumbers));
percent_of_filled = length(N_filled)/length(Recon3D.rxnECNumbers);

%% Get reactions abundances
rxn_abundances = mapAbundanceToReactions(Recon3D, abundance);

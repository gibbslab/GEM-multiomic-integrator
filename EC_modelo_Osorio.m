close all force; clear variables; clc
%% Startup the COBRA Toolbox
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver('glpk', 'all');
%% Load the model
model=readCbModel('Astrocyte_Osorio2019.xml');

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


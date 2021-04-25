close all force; clear variables; clc
%% Startup the COBRA Toolbox
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver('glpk', 'all');
%% Load the model
modelo=readCbModel('Astrocyte_Osorio2019.xml');

% Go to https://www.uniprot.org/uploadlists/ select Entrez to UniprotKB
% Enable the EC numbers column and download as excel
ECnumbers_osorio=readtable('ECnumbers_Osorio.xlsx');
Var_rem=table2array(ECnumbers_osorio(:,'ECNumber'));


for i=1:2405
    Secuencia=strcat(['x(',num2str(i),')'],'');
    modelo.rules =  strrep(modelo.rules, Secuencia,Var_rem(i));
end


%remp = '( )|';
%modelo.rules =  strrep(modelo.rules, remp ,'');
%modelo.rules =  strrep(modelo.rules, '\)','');

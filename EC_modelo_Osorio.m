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
ECNumbers = strrep(ECNumbers,'TCDB:','');
ECNumbers = regexprep(ECNumbers, '\<(?!EC:|^\>)([0-9A-Z]+.(([0-9A-Z]|-)+.)*([0-9A-Z]|-)+)','x($1)');

unique_ECNumbers = {};
sizeEcNumbers = 0;
genes = regexp(model.rules,'(?<=x\()([0-9]+)(?=\))','match');
model.Prules = model.rules;
for i=1:numel(genes)
    indexes = unique(str2double(genes{i}));
    for j=1:numel(indexes)
        index = indexes(j);
        pattern = ['x(' num2str(index) ')'];
        % We need to delete the entry and one of the logic characters around it
        if isempty(ECNumbers{index})
            regExp = regexptranslate('escape', pattern);
            model.Prules(i) = regexprep(model.Prules(i), [regExp ' *[&|] *'], '');
            model.Prules(i) = regexprep(model.Prules(i), ['( *[&|] *)?' regExp], '');
        else
            % Create the table of EC numbers
            match = strcmp(ECNumbers{index}, unique_ECNumbers);
        
            if isempty(match) || sum(match) == 0 
                sizeEcNumbers = sizeEcNumbers + 1;
                ECNumberIndex = sizeEcNumbers;
                % Save in the list of unique ECNumbers
                unique_ECNumbers{sizeEcNumbers} = ECNumbers{index};
            else
                ECNumberIndex = find(match);
            end
            replacement = ['x(' num2str(ECNumberIndex) ')'];
            model.Prules(i) = strrep(model.Prules(i), pattern, replacement);
        end
    end
end

model.Prules = regexprep(model.Prules,'((\||&) *\( *\) *)*$',''); % Remove emty parenthesis at the end
model.Prules = regexprep(model.Prules,' *\( *\) *(\||&)?',''); % Remove emty parenthesis at the start and in the midle
model.ECNumbers = unique_ECNumbers;

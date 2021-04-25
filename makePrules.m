function model = makePrules(model)
%makePrules This function takes model.rxnECNumbers and makes a list with
%the protein rules (model.Prules) and a list of unique EC numbers (model.ECNumbers)
%
% USAGE:
%   model = makePrules(model)
%
% INPUT:
%   model:       cobra model structure
%
% OUTPUT:
%   model:   cobra model with aditional properties Prules and ECNumbers
%
% AUTHORS: Nicolas Mendoza-Mejia, Apr 2020

rxnECNumbers = regexp(model.rxnECNumbers, '\<(?!EC:|^\>)([0-9A-Z]+.(([0-9A-Z]|-)+.)*([0-9A-Z]|-)+)','match');

nRxns = numel(model.rxns);
Prules = cell (nRxns, 1);
ECNumbers = {};
sizeEcNumbers = 0;

for i=1:nRxns

    nECNum = numel(rxnECNumbers{i});
    
    % This is done to avoid adding non-indexed things from rxnECNumbers
    if (nECNum > 0)
        Prules{i} = model.rxnECNumbers{i};
    end
    
    for j=1:nECNum
        match = strcmp(rxnECNumbers{i}{j}, ECNumbers);
        
        if isempty(match) || sum(match) == 0 
            sizeEcNumbers = sizeEcNumbers + 1;
            ECNumberIndex = sizeEcNumbers;
            % Save in the list of unique ECNumbers
            ECNumbers{sizeEcNumbers} = rxnECNumbers{i}{j};
        else
            ECNumberIndex = find(match);
        end
        
        % Modify the protein rules
        pattern = regexptranslate('escape',rxnECNumbers{i}{j});
        replacement = regexptranslate('escape',['x(' num2str(ECNumberIndex) ')']);
        Prules{i} = regexprep(Prules{i}, pattern, replacement, 'once');
    end
   
end

model.Prules = Prules;
model.ECNumbers = ECNumbers;

end


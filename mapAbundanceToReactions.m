function result = mapAbundanceToReactions(modelECNumbers, ProteinECNumbers, ProteinAbundance)
%MAPABUNDANCETOREACTIONS This function maps the protein abundances to 
% the reactions in the model.

%% 3. Find matches between model and abundance
rxnECNumbers = regexp(modelECNumbers, '\<(?!EC:|^\>)([0-9A-Z]+.([0-9A-Z]+.)*([0-9A-Z]|-)+)','match');
rxnECNumbersLogic = regexp(modelECNumbers, '(\||&)','match');

rxn_n = numel(rxnECNumbers);
result = cell (rxn_n, 1);

f = waitbar(0,'Procesing','Name','Progress',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

for i=1:rxn_n
    waitbar(i/rxn_n,f,sprintf('%12.0f',floor(i/rxn_n*100)))
    [Abundance, index] = findAbundances(rxnECNumbers{i}, ProteinECNumbers, ProteinAbundance);

%% 4 and = min,  or = max
    % Now we have the protein abundance
    % We need to check wether the logic is and/or to use min/max
    nLogic = numel(rxnECNumbersLogic{i});
    nEnzyme = numel(rxnECNumbers{i});
    if nEnzyme > 1
        assert(nEnzyme - 1 == nLogic, 'Wrong logic separator for multiple enzymes in model.rxnECNumbers{%i}:%s', i, modelECNumbers{i});

        rxnAbundance = Abundance(1);
        for j=1:nLogic
            if rxnECNumbersLogic{i}{j} == '&'
                rxnAbundance = min(rxnAbundance, Abundance(j + 1));
            elseif rxnECNumbersLogic{i}{j} == '|'
                rxnAbundance = max(rxnAbundance, Abundance(j + 1));
            end
        end

        result{i} =  rxnAbundance;
    elseif nEnzyme == 1 && Abundance ~= 0
        result{i} = Abundance;
    end
    
end
delete(f)
end


function [abundance, index] = findAbundances(rxnECNumbers, ProteinECNumbers, ProteinAbundance)
nEnzymeNumbers = numel(rxnECNumbers);
abundance = zeros(nEnzymeNumbers, 1);
index = zeros(nEnzymeNumbers, 1);

for i=1:nEnzymeNumbers
    rxnECNumber = rxnECNumbers{i};
    %% Find matches
    if ~contains(rxnECNumbers{i}, ".-")
        found = cellfun(@(x) strcmp(x, rxnECNumber), ProteinECNumbers, 'UniformOutput', false);
        match = cellfun(@(c) any(c(:)), found);
    else
        rxnECNumber = strrep(rxnECNumber, '.-', '.');
        found = cellfun(@(x) contains(x, rxnECNumber), ProteinECNumbers, 'UniformOutput', false);
        match = cellfun(@(c) any(c(:)), found);
    end
    %% If there is a match get its info
    if sum(match) > 0
        matchInexes = find(match);
        abundances = str2double(ProteinAbundance(match));
        [abundance(i), localIndex] = max(abundances);
        index(i) = matchInexes(localIndex);
    end
end 

end


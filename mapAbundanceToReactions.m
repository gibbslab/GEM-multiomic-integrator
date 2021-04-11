function result = mapAbundanceToReactions(model, abundance)
%MAPABUNDANCETOREACTIONS This function maps the protein abundances to 
% the reactions in the model.


%% 1. Get the list of ProteinECNumbers
% Maybe we can check if the proteins are already downloaded to not do that
% again
if isfile('ProteinECNumbers.mat')
    load('ProteinECNumbers.mat');
else
    ProteinECNumbers = getProteinCodes(abundance);
    save('ProteinECNumbers.mat', 'ProteinECNumbers');
end

%% 2. Take rxnECNumbers and format correctly 
% .rxnECNumbers ; = & or = | .- = all in the subgroup
model.rxnECNumbers = regexprep(model.rxnECNumbers,";$","");
model.rxnECNumbers = strrep(model.rxnECNumbers,";"," & ");
model.rxnECNumbers = strrep(model.rxnECNumbers,","," & ");
model.rxnECNumbers = strrep(model.rxnECNumbers,"or","|");
% model.rxnECNumbers = strrep(model.rxnECNumbers,".-","all");
model.rxnECNumbers = strrep(model.rxnECNumbers,"EC:","");
model.rxnECNumbers = strrep(model.rxnECNumbers,"TCDB:","");

%% 3. Find matches between model and abundance
rxnECNumbers = regexp(model.rxnECNumbers, '\<(?!EC:|^\>)([0-9A-Z]+.([0-9A-Z]+.)*([0-9A-Z]|-)+)','match');
rxnECNumbersLogic = regexp(model.rxnECNumbers, '(\||&)','match');

rxn_n = numel(rxnECNumbers);
result = cell (rxn_n, 1);

f = waitbar(0,'Procesing','Name','Progress',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

for i=1:rxn_n
    waitbar(i/rxn_n,f,sprintf('%12.0f',floor(i/rxn_n*100)))
    [Abundance, index] = findAbundances(rxnECNumbers{i}, ProteinECNumbers, abundance.Control_NHA1_veh_tech2);
    IDs = abundance.ID(index(index > 1));

%% 4 and = min,  or = max
    % Now we have the protein abundance
    % We need to check wether the logic is and/or to use min/max
    nLogic = numel(rxnECNumbersLogic{i});
    nEnzyme = numel(rxnECNumbers{i});
    if nEnzyme > 1
        assert(nEnzyme - 1 == nLogic, 'Wrong logic separator for multiple enzymes in model.rxnECNumbers{%i}:%s' ,i ,model.rxnECNumbers{i});

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



function ProteinECNumbers = getProteinCodes (abundance)
import matlab.net.*
import matlab.net.http.*

tic;
IDs = abundance.ID;
IDs_length = numel(abundance.ID);
ProteinECNumbers = cell(IDs_length, 1);

f = waitbar(0,'Please wait...');

parfor i=1:IDs_length
    UniProtURL = ['https://www.uniprot.org/uniprot/' IDs{i} '.txt'];
    txt = webread(UniProtURL); 
    
    % Get the Enzyme Commission numbers
    ECNumbers = regexp(txt, 'DE[ ]*EC=(?<ECNumbers>([0-9]+.([0-9]+.)*([0-9]|-)+))','names');
    % Get the Transporter Classification numbers 
    TCDBNumbers = regexp(txt, 'DR[ ]*TCDB; (?<TCDBNumbers>([0-9A-Z]+.([0-9A-Z]+.)*([0-9A-Z]|-)+))','names');
    

   ProteinECNumbers{i} = [
       reshape(struct2cell(ECNumbers), numel(ECNumbers), 1)
       reshape(struct2cell(TCDBNumbers), numel(TCDBNumbers), 1)
       ];
    
    waitbar(i/IDs_length, f,['Downloading (' num2str(i) '/' num2str(IDs_length) ')...']);
end
close(f)
toc;


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

